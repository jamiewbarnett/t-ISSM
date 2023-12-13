#Copy files to specific folders

now=$(date +"%Y_%m_%d")



rm *.tar.gz; rm *.bin; rm *.asv; rm *.errlog; rm *.outlog; rm *.queue; rm *.toolkits; 


git add --all -- ':!Model_Data/*' ':!Outputs/*' ':!.DS_Store'
git add Model_Data/GIS/*
git add Outputs/outputs.txt
git commit -a -m "Lazy update: ${now}"
git pull
git push

