#Copy files to specific folders

now=$(date +"%Y_%m_%d")



rm *.tar.gz; rm *.bin; 

git add --all -- ':!Model_Data/*' ':!Models/*'
git commit -a -m "Lazy update: ${now}"
git push

