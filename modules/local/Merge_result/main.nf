process Merge_result { 
  publishDir params.output+'/merge/', mode: 'copy'
  
  input:
    file bam
  
  output:
    file "finale.txt", emit: Merge_result
  
  shell:
  '''
  ## -- Differancial expression analyse ------------------------------------------------- ##
  files=(*)
  awk '{print $1}' ${files[0]} > AAAA.txt

  for file in *_intersect.txt; do
    awk '{print $2}' $file > ${file}.tmp
    rm $file
    mv ${file}.tmp $file 
    tail -n +2 "$file" > "$file.tmp" && mv "$file.tmp" "$file"
    echo "${file%%.*}" > ${file}.name
    cat ${file}.name ${file} > ${file}.tmp && mv ${file}.tmp ${file}
    rm ${file}.name
  done

  paste -d "\t" * > finale.txt
  rm AAAA.txt
  '''
}
