# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2019-2025, The OpenROAD Authors

add_library(bookshelf
    bookshelfout_impl.cpp
    bookshelfout.cpp
)

target_include_directories(bookshelf
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
        ${TCL_INCLUDE_PATH}
)
target_link_libraries(bookshelf
    db
    utl_lib
)

set_target_properties(bookshelf
  PROPERTIES
    # python requirement
    POSITION_INDEPENDENT_CODE ON
)

