Errors with setting working directory
There can be issues in getting the right path, this seems to often fix the issue
1. go to the python "site-packages" folder
2. create a new text file
3. add the absolute path to "Flux_Code"
4. save with the .pth extension
Alternatively there can be issues if your cwd starts at Flux_Code, adding 'if not Path(os.getcwd()).parts[-1] == "Flux_Code":' right above 'path_to_FC = ""' can fix this