import pandas as pd
import cobra
from cobra.io import read_sbml_model
from cobra.flux_analysis import pfba

# Load the COBRA model
model = cobra.io.read_sbml_model("iTP252_irreversible_model.xml")

# Set the bounds of specific reactions to 0
reactions_bound_0 = [
    "EX_cpd00020_e0_b", "rxn01512_c0_b", "rxn01513_c0_b", "rxn01127_c0_b",
    "rxn00412_c0_b", "rxn00410_c0_b", "rxn08192_c0_b", "rxn05148_c0_b",
    "rxn00119_c0_b", "rxn00770_c0_b", "rxn01517_c0_b", "rxn00225_c0_f",
    "rxn00097_c0_b", "rxn00392_c0_b", "rxn02314_c0_b", "rxn01100_c0_f",
    "rxn00216_c0_b", "rxn00077_c0_b", "rxn00364_c0_b", "rxn01673_c0_b",
    "rxn01219_c0_b", "rxn00237_c0_b", "rxn01678_c0_b", "rxn00515_c0_b",
    "rxn01353_c0_b", "rxn02155_c0_b", "rxn00409_c0_b", "rxn02517_c0_b",
    "rxn00117_c0_b", "rxn00839_c0_b", "rxn00190_c0_b", "EX_cpd00027_e0_f",
    "EX_cpd00023_e0_f", "EX_cpd00129_e0_f", "EX_cpd00053_e0_f", "EX_cpd00184_e0_f"
]

for reaction_id in reactions_bound_0:
    reaction = model.reactions.get_by_id(reaction_id)
    reaction.lower_bound = 0
    reaction.upper_bound = 0

# Set specific bounds for another reaction
model.reactions.get_by_id("EX_cpd00027_e0_b").upper_bound = 0.75
model.reactions.get_by_id("EX_cpd00027_e0_b").lower_bound = 0.75

# Read reactions and protein costs from the Excel file
reaction_data = pd.read_excel("kcat_mw.xlsx", usecols=[0, 1], names=["Reactions", "Protein_cost"])

# Convert DataFrame to dictionary for easy lookup
reaction_protein_costs = dict(zip(reaction_data["Reactions"], reaction_data["Protein_cost"]))

# Get the reactions from the model and calculate the sum of fluxes multiplied by protein costs
total_protein_cost_expression = sum(
    model.reactions.get_by_id(reaction).flux_expression * reaction_protein_costs[reaction]
    for reaction in reaction_protein_costs.keys() if reaction in model.reactions
)

# Set a constraint on the total protein content to 700 (if you want to fix this as a constraint)
total_protein_cost_constraint = model.problem.Constraint(total_protein_cost_expression, lb=0, ub=90000)
model.add_cons_vars([total_protein_cost_constraint])

# Fix the flux of bio1_biomass to a specific value, e.g., 0.0231
biomass_reaction = model.reactions.get_by_id('bio1_biomass')
biomass_reaction.lower_bound = 0.0231
biomass_reaction.upper_bound = 0.0231

# Set the objective to maximize the total protein cost
model.objective = model.problem.Objective(total_protein_cost_expression, direction='min')

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

# Print the total protein cost
print("Total protein cost:", total_protein_cost)

# Write the DataFrame to an Excel file, noting that this is the pFBA solution
output_filename = "find_lowest_protein.xlsx"
with pd.ExcelWriter(output_filename, engine='openpyxl') as writer:
    reactions_df.to_excel(writer, sheet_name='Reactions Flux', index=False)
    atp_fluxes_df.to_excel(writer, sheet_name='ATP Producing Reactions', index=False)

print(f"Output written to {output_filename}")

# Print the biomass flux
print("Biomass flux:", biomass_flux)
