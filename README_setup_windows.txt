# GitHub Setup Guide (Windows)

## Step 1: Install Git

Download Git for Windows from:
https://git-scm.com/

Install it with default settings.

## Step 2: Generate a GitHub Token

1. Go to: https://github.com/settings/personal-access-tokens
2. Click **"Generate new token (fine-grained)"**
3. Fill in:
   - Repository access: Select only `ProjetBioinfo`
   - Permissions: Enable `Read and Write` under "Repository permissions"
   - Set expiration to 30 days or as needed
4. Click "Generate token"
5. **Copy the token and save it** (you will not see it again)

## Step 3: Clone the Repository

Open Git Bash and run:

git clone https://github.com/USERNAME/ProjetBioinfo.git

When prompted:
- Enter your GitHub username
- Paste your token as the password

## Step 4: Work with Git

cd ProjetBioinfo

To add and push changes:

git add .
git commit -m "Your message"
git push

To pull updates:

git pull

## Notes

- Use `git config --global credential.helper wincred` to cache credentials.
- Never share your token. You can revoke it in your GitHub settings.
