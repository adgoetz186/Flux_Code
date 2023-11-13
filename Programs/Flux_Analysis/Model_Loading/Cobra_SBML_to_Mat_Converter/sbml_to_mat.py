import cobra

model = cobra.io.read_sbml_model("/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/Models/sbml_models/MODEL1603150001_url.xml")
cobra.io.save_matlab_model(model, "/Users/agoetz/PycharmProjects/pythonProject2/Flux_Code/Data/Models/mat_models/recon2_2.mat")