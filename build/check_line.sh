#!/bin/bash

# 遍历所有 good_query_sequences_ 开头的文件
for good_query_file in good_query_sequences_*.fasta; do
    # 提取后缀
    suffix="${good_query_file#good_query_sequences_}"

    # 对应的 ref_sequences_ 文件
    ref_file="ref_sequences_$suffix"

    # 检查 ref_sequences_ 文件是否存在
    if [[ -f "$ref_file" ]]; then
        # 获取两个文件的行数
        good_query_lines=$(wc -l < "$good_query_file")
        ref_lines=$(wc -l < "$ref_file")

        # 比较行数
        if [ "$good_query_lines" -eq "$ref_lines" ]; then
            echo "文件 $good_query_file 和 $ref_file 行数相同。"
        else
            echo "文件 $good_query_file 和 $ref_file 行数不同。"
        fi
    else
        echo "没有找到与 $good_query_file 匹配的 $ref_file。"
    fi
done

