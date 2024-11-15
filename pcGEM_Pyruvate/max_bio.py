import pandas as pd
import cobra
from cobra.io import read_sbml_model
from cobra.flux_analysis import pfba

# Load the COBRA model
model = cobra.io.read_sbml_model("iTP251_irreversible_model.xml")

reactions_bound_0 = [
    "EX_cpd00027_e0_b", "rxn01512_c0_b", "rxn01513_c0_b", "rxn01127_c0_b",
    "rxn00412_c0_b", "rxn00410_c0_b", "rxn08192_c0_b", "rxn05148_c0_b",
    "rxn00119_c0_b", "rxn00770_c0_b", "rxn01517_c0_b", "rxn00225_c0_f",
    "rxn00097_c0_b", "rxn00392_c0_b", "rxn02314_c0_b",
    "rxn00077_c0_b", "rxn00364_c0_b", "rxn01673_c0_b",
    "rxn01219_c0_b", "rxn00237_c0_b", "rxn01678_c0_b", "rxn00515_c0_b",
    "rxn01353_c0_b", "rxn02155_c0_b", "rxn00409_c0_b", "rxn02517_c0_b",
    "rxn00117_c0_b", "rxn00839_c0_b", "rxn00190_c0_b", "EX_cpd00020_e0_f",
    "EX_cpd00138_e0_b", "EX_cpd00138_e0_f"
]

for reaction_id in reactions_bound_0:
    reaction = model.reactions.get_by_id(reaction_id)
    reaction.lower_bound = 0
    reaction.upper_bound = 0

aa_reactions = [
    "EX_cpd00041_e0_b", "EX_cpd00023_e0_b", "EX_cpd00065_e0_b", "EX_cpd00084_e0_b",
    "EX_cpd00053_e0_b", "EX_cpd00132_e0_b", "EX_cpd00107_e0_b", "EX_cpd00129_e0_b",
    "EX_cpd00322_e0_b", "EX_cpd00039_e0_b", "EX_cpd00054_e0_b", "EX_cpd00033_e0_b",
    "EX_cpd00035_e0_b", "EX_cpd00156_e0_b", "EX_cpd00161_e0_b", "EX_cpd00066_e0_b",
    "EX_cpd00053_e0_b", "EX_cpd00041_e0_b", "EX_cpd00069_e0_b", "EX_cpd00051_e0_b",
    "EX_cpd00060_e0_b", "EX_cpd00119_e0_b", "EX_cpd00159_e0_b", "EX_cpd00221_e0_b"
    
]
for reaction_id in aa_reactions:
    reaction = model.reactions.get_by_id(reaction_id)
    reaction.lower_bound = 0
    reaction.upper_bound = 0.3


model.reactions.get_by_id("EX_cpd00020_e0_b").lower_bound = 0.3
model.reactions.get_by_id("EX_cpd00020_e0_b").upper_bound = 0.3

# Fix the flux of bio1_biomass to a specific value, e.g., 0.0231 (if necessary)
biomass_reaction = model.reactions.get_by_id('bio1_biomass')

# Set the objective to maximize the biomass reaction flux
model.objective = biomass_reaction

# Run parsimonious FBA (pFBA) to maximize biomass
pfba_solution = pfba(model)

# Extract the flux value for the biomass reaction from the solution
biomass_flux = pfba_solution.fluxes.get('bio1_biomass')

# Extract the flux values for all reactions using the pFBA solution
fluxes = pfba_solution.fluxes

# Calculate the total protein cost using flux values from the pFBA solution (optional)
reaction_data = pd.read_excel("kcat_mw.xlsx", usecols=[0, 1], names=["Reactions", "Protein_cost"])
reaction_protein_costs = dict(zip(reaction_data["Reactions"], reaction_data["Protein_cost"]))

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
output_filename = "find_maximum_biomass.xlsx"
with pd.ExcelWriter(output_filename, engine='openpyxl') as writer:
    reactions_df.to_excel(writer, sheet_name='Reactions Flux', index=False)
    atp_fluxes_df.to_excel(writer, sheet_name='ATP Producing Reactions', index=False)

print(f"Output written to {output_filename}")

# Print the biomass flux
print("Biomass flux:", biomass_flux)
