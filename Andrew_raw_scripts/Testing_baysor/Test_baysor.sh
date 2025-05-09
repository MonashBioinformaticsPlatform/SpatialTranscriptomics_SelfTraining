#!/bin/bash
LD_LIBRARY_PATH=""

toml="/homevol/apattison/Code/Spatial/Paper_code_v2/Baysor/CosMx.toml"
out="/oldvol/apattison/Data/Spatial/Comparing_technologies/CosMx_Abud_Baysor/baysor_out/"
mols="/oldvol/apattison/Data/Spatial/Comparing_technologies/CosMx_Abud_Baysor/Run5654_Tumor_A_tx_file_fovs_5-8_test.csv"
cd /pvol/andrew/bin/baysor/bin/baysor/bin

./baysor run -c "${toml}" -o "${out}" "${mols}" ":cell_ID" > "${out}/baysor_fovs_5-8_test.log" 2>&1

