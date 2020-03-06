
# install txt files
set( ${product}_file_dir "${product}/${version}/files" )
file(GLOB txt_file *.txt )
install( FILES ${txt_file} DESTINATION ${${product}_file_dir} )

