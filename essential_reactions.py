import cobra
import pandas as pd

# Load the model
model = cobra.io.read_sbml_model('iTP251.xml')

# Initialize lists to store results
reactions_list = []
biomass_reduction_list = []
pathway_list = []

# Get the original biomass reaction objective value
original_solution = model.optimize()
original_biomass = original_solution.objective_value

# Iterate through each reaction in the model
for reaction in model.reactions:
    # Temporarily knock out the reaction
    with model:
        reaction.knock_out()
        solution = model.optimize()
        
        # Calculate biomass reduction
        biomass_reduction = (original_biomass - solution.objective_value) / original_biomass * 100
        
        # Append data to lists
        reactions_list.append(reaction.id)
        biomass_reduction_list.append(biomass_reduction)
        pathway_list.append(reaction.name)  # Use the reaction name as the pathway

# Create a DataFrame to store the results
results_df = pd.DataFrame({
    'Reaction': reactions_list,
    'Biomass Reduction (%)': biomass_reduction_list,
    'Enzyme': pathway_list
})

# Save the DataFrame to an Excel file
results_df.to_excel('essential_reactions2.xlsx', index=False)

print("Analysis complete. All reactions saved to reactions_analysis_with_names.xlsx")