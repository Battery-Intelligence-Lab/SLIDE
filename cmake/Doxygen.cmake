function(enable_doxygen input output)
  option(ENABLE_DOXYGEN "Enable doxygen doc builds of source" ON)
  if(ENABLE_DOXYGEN)

    # set(DOXYGEN_CALLER_GRAPH YES)
    # set(DOXYGEN_CALL_GRAPH YES)
    # set(DOXYGEN_EXTRACT_ALL YES)

    #set(DOXYGEN_ALPHABETICAL_INDEX NO)
    #set(DOXYGEN_SHOW_NAMESPACES NO)
    #set(DOXYGEN_CREATE_SUBDIRS YES) -> Uses random file names. 

    set(DOXYGEN_GENERATE_TREEVIEW YES)
    set(DOXYGEN_GENERATE_HTML YES)
    set(DOXYGEN_DOT_IMAGE_FORMAT svg)
    set(DOXYGEN_HTML_OUTPUT
      ${PROJECT_BINARY_DIR}/${output})

    # set(DOXYGEN_HTML_EXTRA_STYLESHEET
    #   ${PROJECT_SOURCE_DIR}/../doxygen-awesome-css/doxygenawesome.css)

    find_package(Doxygen REQUIRED dot)
    doxygen_add_docs(doxygen-docs ${PROJECT_SOURCE_DIR}/${input})

  endif()
endfunction()
