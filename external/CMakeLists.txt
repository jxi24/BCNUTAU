CPMFindPackage(
    NAME fmt
    GIT_TAG 7.1.3
    GITHUB_REPOSITORY fmtlib/fmt
)
add_library(fmt::fmt ALIAS fmt)

CPMFindPackage(
    NAME spdlog
    GIT_TAG "v1.8.5"
    GITHUB_REPOSITORY gabime/spdlog
    OPTIONS
        "SPDLOG_COMPILED_LIB ON"
        "SPDLOG_FMT_EXTERNAL ON"
        "SPDLOG_BUILD_SHARED ON"
)
add_library(spdlog::spdlog ALIAS spdlog)

CPMAddPackage(
    NAME Chili
    # VERSION 1.0.1
    # URL https://gitlab.com/spice-mc/Chili/-/archive/v1.0/Chili-v1.0.1.tar.gz
    URL https://gitlab.com/spice-mc/Chili/-/archive/subproject/Chili-main.tar.gz
)
add_library(Chili::Chili ALIAS Chili)
    

# Install testing framework
if(CHILI_ENABLE_TESTING)
    # Catch2
    CPMAddPackage(
        NAME catch2
        GITHUB_REPOSITORY catchorg/Catch2
        GIT_TAG "v2.13.7"
    )
    add_library(Catch2::Catch2 ALIAS Catch2)

    # Trompeloeil
    CPMAddPackage(
        NAME trompeloeil
        GITHUB_REPOSITORY rollbear/trompeloeil
        VERSION 37
    )
endif()
