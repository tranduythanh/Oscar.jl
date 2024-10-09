#!/bin/sh

# Fetch the latest changes from the remote
git fetch origin

# Check if there are any changes to pull
if git status -uno | grep -q 'Your branch is behind'; then
    echo "Changes detected on remote. Pulling changes..."
    
    # Try to pull changes
    if ! git pull --no-ff; then
        echo "Pull failed. Please resolve conflicts manually and then run the script again."
        exit 1
    fi
else
    echo "Local branch is up to date."
fi

# Show current status
git status

# Add all changes
git add .

# Show status after adding changes
git status

# Commit changes
git commit -m "add new functionality"

# Show status after commit
git status

# Push changes to remote
git push origin feature/mr-hiep

echo "done ---------<>"
