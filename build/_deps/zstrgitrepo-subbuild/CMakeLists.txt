# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.16.3)

# We name the project and the target for the ExternalProject_Add() call
# to something that will highlight to the user what we are working on if
# something goes wrong and an error message is produced.

project(zstrgitrepo-populate NONE)

include(ExternalProject)
ExternalProject_Add(zstrgitrepo-populate
                     "UPDATE_DISCONNECTED" "False" "GIT_REPOSITORY" "https://github.com/mateidavid/zstr" "GIT_TAG" "755da7890ea22478a702e3139092e6c964fab1f5"
                    SOURCE_DIR          "/home/user_home/ylf/someGit/strobealign-0.12.0-timer/build/_deps/zstrgitrepo-src"
                    BINARY_DIR          "/home/user_home/ylf/someGit/strobealign-0.12.0-timer/build/_deps/zstrgitrepo-build"
                    CONFIGURE_COMMAND   ""
                    BUILD_COMMAND       ""
                    INSTALL_COMMAND     ""
                    TEST_COMMAND        ""
                    USES_TERMINAL_DOWNLOAD  YES
                    USES_TERMINAL_UPDATE    YES
)
