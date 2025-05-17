# GitHub Setup Guide

## Step 1: Install Git

Download Git for Windows from:
https://git-scm.com/

Install it with default settings.

## Step 2: Download gitHub desktop

https://desktop.github.com/download/

Login github using your account that has been granted access to the GitHub
Choose the Repository and open it in VS code (option should be visible in the right pannel)

## Step 3: Clone the Repository

choose the "clone repository" option on first opening, it will create a folder in the GitHub files.

## Step 4: Work with Git

You can eddit things in VS normally. Once you're happy with the changes, save them, 
go to the GitHub app and add a short commentary about your changes. Hit 'Commit to main', 
then hit the 'Push' button on the top of the window. Your changes have now been uploaded to github.

To add and push changes manually in the terminal:

git add .
git commit -m "Your message"
git push

To pull updates:

git pull

## Notes

- Use `git config --global credential.helper wincred` to cache credentials.
- Never share your token. You can revoke it in your GitHub settings.

---

## Step 5: Handling Large Data Files (Google Drive)

GitHub does not allow pushing files over 100MB. Use Google Drive for large datasets.

### Recommended Setup

1. Create a folder in Google Drive (e.g., `ProjetBioinfo_data`)
2. Share it with collaborators, or add a shortcut from the shared file in your main Drive
3. Sync it locally using the Google Drive desktop app
4. Create a shortcut named 'data' in the cloned repository on your machine that links to the data folder on your drive.

Make sure to add `data/` to your `.gitignore` (it should already be in there):

```
data/
```

This prevents Git from trying to track large files.

### Summary

- Code is tracked via GitHub
- Data is stored in Drive and manually managed per collaborator
