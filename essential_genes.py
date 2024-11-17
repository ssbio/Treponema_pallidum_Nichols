import cobra
import numpy as np

# Load the metabolic model from the current directory
model = cobra.io.read_sbml_model("iTP251.xml")

# Define the objective function to maximize growth rate
model.objective = model.reactions.get_by_id("bio1_biomass")

# Perform single gene deletions and simulate growth
essential_genes = []
for gene in model.genes:
    # Create a copy of the model with the gene deleted
    model_copy = model.copy()
    model_copy.genes.get_by_id(gene.id).knock_out()

    # Simulate growth and check if it is reduced by 90% or more
    wt_sol = model.optimize()
    del_sol = model_copy.optimize()
    if del_sol.objective_value / wt_sol.objective_value <= 0.10:
        essential_genes.append(gene.id)

# Write the list of essential genes to a text file
with open("essential_genes.txt", "w") as f:
    for gene in essential_genes:
        f.write(gene + "\n")

# Print a message indicating where the output was written
print("Essential genes written to essential_genes.txt.")