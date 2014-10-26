function build_nova-daq-offline() {
  prep_build nova-daq-offline ${1} ${2}:${build_type} || return 0
  ./build_nova-daq-offline.sh ${product_topdir} ${2} ${build_type} ${maketar} >& "${logfile}"
}

# Local Variables:
# mode: sh
# eval: (sh-set-shell "bash")
# End:
