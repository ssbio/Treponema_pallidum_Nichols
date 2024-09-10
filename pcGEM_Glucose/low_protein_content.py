import pandas as pd
import cobra
from cobra.io import read_sbml_model
from cobra.flux_analysis import pfba

# Load the COBRA model
model = cobra.io.read_sbml_model("iTP252_irreversible_model.xml")

# Set the bounds of specific reactions to 0
reactions_bound_0 = [
    "EX_cpd00020_e0_b",
    "EX_cpd00023_e0_f", "EX_cpd00129_e0_f", "EX_cpd00053_e0_f", "EX_cpd00184_e0_f",
]

for reaction_id in reactions_bound_0:
    reaction = model.reactions.get_by_id(reaction_id)
    reaction.lower_bound = 0
    reaction.upper_bound = 0

# Set specific bounds for another reaction
model.reactions.get_by_id("EX_cpd00027_e0_r").upper_bound = 0.75
model.reactions.get_by_id("EX_cpd00027_e0_r").lower_bound = 0.75

# Read reactions and protein costs from the Excel file
reaction_data = pd.read_excel("kcat_mw_new.xlsx", usecols=[0, 1], names=["Reactions", "Protein_cost"])

# Convert DataFrame to dictionary for easy lookup
reaction_protein_costs = dict(zip(reaction_data["Reactions"], reaction_data["Protein_cost"]))

# Get the reactions from the model and calculate the sum of fluxes multiplied by protein costs
total_protein_cost_expression = sum(
    model.reactions.get_by_id(reaction).flux_expression * reaction_protein_costs[reaction]
    for reaction in reaction_protein_costs.keys() if reaction in model.reactions
)

# Fix the protein cost constraint to 77
total_protein_cost_constraint = model.problem.Constraint(total_protein_cost_expression, lb=96.62, ub=96.62)
model.add_cons_vars([total_protein_cost_constraint])

# Fix the flux of bio1_biomass to 0.0231 for the constraint
biomass_reaction = model.reactions.get_by_id('bio1_biomass')

# Set the objective to maximize biomass
model.objective = biomass_reaction

# Run parsimonious FBA (pFBA)
pfba_solution = pfba(model)

# Extract the flux value for the biomass reaction from the solution
biomass_flux = pfba_solution.fluxes.get('bio1_biomass')

# Extract the flux values for all reactions using the pFBA solution
fluxes = pfba_solution.fluxes

# Calculate the total protein cost using flux values from the pFBA solution
total_protein_cost = sum(
    fluxes[reaction] * reaction_protein_costs.get(reaction, 0)
    for reaction in fluxes.index if reaction in reaction_protein_costs
)

# Create a DataFrame for the reaction fluxes
reactions_df = pd.DataFrame({
    "Reaction ID": fluxes.index,
    "Flux": fluxes.values
})

# Find all reactions that produce ATP (cpd00002[c0])
atp_metabolite = model.metabolites.get_by_id('cpd00002[c0]')
atp_producing_reactions = [reaction for reaction in atp_metabolite.reactions if atp_metabolite in reaction.products]

# Create a DataFrame for ATP-producing reactions and their fluxes
atp_fluxes_df = pd.DataFrame({
    "Reaction ID": [reaction.id for reaction in atp_producing_reactions],
    "Flux": [fluxes[reaction.id] for reaction in atp_producing_reactions]
})

# Find all reactions that regenerate NAD (cpd00003[c0])
nad_metabolite = model.metabolites.get_by_id('cpd00003[c0]')
nad_regenerating_reactions = [reaction for reaction in nad_metabolite.reactions if nad_metabolite in reaction.products]

# Create a DataFrame for NAD-regenerating reactions and their fluxes
nad_fluxes_df = pd.DataFrame({
    "Reaction ID": [reaction.id for reaction in nad_regenerating_reactions],
    "Flux": [fluxes[reaction.id] for reaction in nad_regenerating_reactions]
})


# Print the total protein cost
print("Total protein cost:", total_protein_cost)

# Write the DataFrame to an Excel file, noting that this is the pFBA solution
output_filename = "lowest_protein_flux_distribution.xlsx"
with pd.ExcelWriter(output_filename, engine='openpyxl') as writer:
    reactions_df.to_excel(writer, sheet_name='Reactions Flux', index=False)
    atp_fluxes_df.to_excel(writer, sheet_name='ATP Producing Reactions', index=False)
    nad_fluxes_df.to_excel(writer, sheet_name='NAD Regenerating Reactions', index=False)

print(f"Output written to {output_filename}")

# Print the biomass flux
print("Biomass flux:", biomass_flux)
