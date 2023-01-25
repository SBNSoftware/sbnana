files=`samweb -e sbn list-definition-files $1`
for f in $files
do
  samweb -e sbn get-file-access-url --schema=xroot $f
done
