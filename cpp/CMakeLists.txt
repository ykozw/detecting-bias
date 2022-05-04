cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
include(FetchContent)

# ソリューションエクスプローラーにcmake用のフォルダを追加する
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER "CMake")
set(CMAKE_CXX_STANDARD 23)

# ソースファイル
set(SRC_FILES
    welch.cpp
    CMakeLists.txt
    )

# プロジェクト作成
project(biasdetector VERSION 1.0)
add_executable(biasdetector ${SRC_FILES})
install (TARGETS biasdetector DESTINATION ${CMAKE_SOURCE_DIR}/bin)

# HACK: ビルドと同時にinstallする
add_custom_command(
  TARGET biasdetector 
    POST_BUILD COMMAND ${CMAKE_COMMAND} --install ${CMAKE_BINARY_DIR} --config $<CONFIG>
)