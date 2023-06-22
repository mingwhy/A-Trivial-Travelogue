@files=glob('./*tab.gz');
#print @files;
for(@files){
    `gzip -d $_`;
}

