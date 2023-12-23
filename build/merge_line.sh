#!/bin/bash

# 遍历所有 good_query_sequences_ 开头的文件
for good_query_file in good_query_sequences_*.fasta; do
    # 提取后缀
    suffix="${good_query_file#good_query_sequences_}"

    # 对应的 ref_sequences_ 文件
    ref_file="ref_sequences_$suffix"

    cat $good_query_file >> merge_query.fasta
    cat $ref_file >> merge_ref.fasta

done

