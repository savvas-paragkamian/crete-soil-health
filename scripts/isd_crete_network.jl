using FlashWeave

data_path = string("results/")

files = readdir(data_path)
for file in files
    file_path = string(data_path, '/', file)
    print(file_path)
#    network = learn_network(file_path, sensitive=true, heterogeneous=false, n_obs_min=10)
#    saveloc = string(data_path, '/', file[begin:end-5], ".gml")
#    save_network(saveloc, network)
end
