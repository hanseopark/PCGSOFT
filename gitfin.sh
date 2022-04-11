DATE=${date}

git add .
git commit -m "changes made on $DATE"

git push origin master

oascript -e 'display noficiation "pushed to remote' with title 'SUCCESS"'
