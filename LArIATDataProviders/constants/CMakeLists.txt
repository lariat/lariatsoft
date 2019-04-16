
# install txt files
set( ${product}_file_dir "${product}/${version}/files" )
file(GLOB txt_file *.txt )
install( FILES ${txt_file} DESTINATION ${${product}_file_dir} )

# also install ROOT files
file(GLOB root_file *.root )
install( FILES ${root_file} DESTINATION ${${product}_file_dir} )


