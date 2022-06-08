project (bam_mergeCore CXX)

add_executable(bam_merge main.cpp)

target_link_libraries(bam_merge zlib bam stdc++fs)
