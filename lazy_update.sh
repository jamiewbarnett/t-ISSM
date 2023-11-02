#Copy files to specific folders

now=$(date +"%Y_%m_%d")



rm *.tar.gz; rm *.bin; rm *.asv 

git pull

git add --all -- ':!Model_Data/*' ':!Outputs/*' ':!.DS_Store'
git add Model_Data/GIS/*
git commit -a -m "Lazy update: ${now}"
git push

