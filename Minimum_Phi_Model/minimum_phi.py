import cobra
import gurobipy as gp
from gurobipy import GRB
import pandas as pd
import matplotlib.pyplot as plt


def read_excel_data(xlsx_path):
    xls = pd.ExcelFile(xlsx_path)
    kcat_df = pd.read_excel(xls, 'Kcat', index_col=0)
    mw_df = pd.read_excel(xls, 'MW', index_col=0)
    return kcat_df, mw_df


def update_reaction_bounds(model, reactions_bound_10):
    """Updates reaction bounds in the model based on specified lists."""
    for rxn_id in reactions_bound_10:
        model.reactions.get_by_id(rxn_id).bounds = (0, 10)
    model.reactions.get_by_id("bio1_biomass").bounds = (0.73338, 0.73338)

def optimize_phi(model, pi_values):
    m = gp.Model("minimize_phi")
    fluxes = {rxn.id: m.addVar(lb=0, ub=1000, name=rxn.id) for rxn in model.reactions}
    m.update()
    
    # Add mass balance constraints for each metabolite
    for met in model.metabolites:
        m.addConstr(
            gp.quicksum(fluxes[rxn.id] * rxn.metabolites[met] for rxn in met.reactions if rxn.id in fluxes) == 0, 
            name=f"mass_balance_{met.id}"
        )
    # Add a constraint that forces non-zero flux through an essential reaction
    m.addConstr(fluxes["bio1_biomass"] >= 0.73338, "demand_biomass")

    # Define the objective to minimize Φ with a small weight on the sum of fluxes to minimize them as a secondary objective
    primary_objective = gp.quicksum((fluxes[rxn_id] / 1000) * pi for rxn_id, pi in pi_values.items())
    secondary_objective = gp.quicksum((fluxes[rxn_id] / 1000) for rxn_id in fluxes) * 0.001  # Adjust the weight as necessary
    m.setObjective(primary_objective + secondary_objective, GRB.MINIMIZE)

    m.optimize()
    
    if m.status == GRB.OPTIMAL:
        return m.objVal, {rxn: fluxes[rxn].x for rxn in fluxes}
    else:
        return float('inf'), {}

# Load the model
model = cobra.io.read_sbml_model("iTP251_irreversible_model.xml")

# Update model bounds

reactions_bound_10 = [
    "EX_cpd00107_e0_b", "EX_cpd00117_e0_b", "EX_cpd00039_e0_b", "EX_cpd00276_e0_b",
    "EX_cpd00041_e0_b", "EX_cpd00069_e0_b", "EX_cpd00156_e0_b", "EX_cpd00023_e0_b",
    "EX_cpd00161_e0_b", "EX_cpd00027_e0_b", "EX_cpd03847_e0_b", "EX_cpd00065_e0_b",
    "EX_cpd00060_e0_b", "EX_cpd00054_e0_b", "EX_cpd00322_e0_b", "EX_cpd00051_e0_b",
    "EX_cpd00393_e0_b", "EX_cpd00132_e0_b", "EX_cpd00119_e0_b", "EX_cpd00367_e0_b",
    "EX_cpd00035_e0_b", "EX_cpd00084_e0_b", "EX_cpd00066_e0_b", "EX_cpd00129_e0_b",
    "EX_cpd00075_e0_b", "EX_cpd00007_e0_b", "EX_cpd00020_e0_b", "EX_cpd00221_e0_b",
    "EX_cpd00138_e0_b"
]

update_reaction_bounds(model, reactions_bound_10)


kcat_df, mw_df = read_excel_data("Kcat_MW_1000simulation_input.xlsx")  # Update this path

phi_results = []
all_fluxes = {}
min_phi = float('inf')
optimal_flux_values = None
optimal_simulation = None

# Perform optimizations for each simulation set
for sim_num in range(1, 101):
    sim_str = f'Simulation {sim_num}:'  # Assuming the colon is part of the header based on the error
    if sim_str in kcat_df.columns and sim_str in mw_df.columns:
        sim_kcat = kcat_df[sim_str].dropna()
        sim_mw = mw_df[sim_str].dropna()
        pi_values = {rxn: sim_mw[rxn] / sim_kcat[rxn] for rxn in sim_kcat.index.intersection(sim_mw.index)}
        
        phi, fluxes = optimize_phi(model, pi_values)

    
    phi_results.append(phi)
    all_fluxes[sim_num] = fluxes
    
    if phi < min_phi:
        min_phi = phi
        optimal_flux_values = fluxes
        optimal_simulation = sim_num
# Write the results to an Excel file
with pd.ExcelWriter('optimization_results.xlsx', engine='openpyxl') as writer:  # Update this path
    pd.DataFrame({'Simulation': range(1, 101), 'Phi Value': phi_results}).to_excel(writer, sheet_name='Phi Values')
    if optimal_flux_values:
        pd.DataFrame({'Reaction ID': list(optimal_flux_values.keys()), 'Flux': list(optimal_flux_values.values())}).to_excel(writer, sheet_name=f'Optimal Fluxes Sim {optimal_simulation}')

print("Optimization completed. Results saved.")

# Now, create a histogram of the Phi values
plt.figure(figsize=(10, 6))
plt.hist(phi_results, bins=20, color='skyblue', edgecolor='black')
plt.title('Distribution of Phi Values Across 100 Simulations')
plt.xlabel('Phi Value')
plt.ylabel('Frequency')
plt.grid(axis='y', alpha=0.75)

# Save the histogram to a file
plt.savefig('phi_values_histogram.png', dpi=300)

# Optionally, show the histogram in a window (this line can be omitted if running in a non-interactive environment)
plt.show()

print("Optimization completed. Results and histogram saved.")