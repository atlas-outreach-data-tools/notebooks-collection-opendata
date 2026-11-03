#!/bin/bash

ENV_FILE=$1

if [ -z "$ENV_FILE" ]; then
    echo "Usage: ./create_environment.sh environment.yml"
    exit 1
fi

# read environment name and Python version

ENV_NAME=$(grep '^name:' "$ENV_FILE" | awk '{print $2}')
PYTHON_VERSION=$(grep '^python:' "$ENV_FILE" | awk '{print $2}')

echo "Environment: $ENV_NAME"
echo "Python:      $PYTHON_VERSION"

# check that we are not currently in the environment

CURRENT_ENV=$(pyenv version-name)

if [ "$CURRENT_ENV" = "$ENV_NAME" ]; then
    echo ""
    echo "ERROR: You are currently inside '$ENV_NAME'."
    echo "Please deactivate it first:"
    echo ""
    echo "    pyenv deactivate"
    echo ""
    exit 1
fi

# remove existing environment

if pyenv versions --bare | grep -q "^${ENV_NAME}$"; then
    echo ""
    echo "ERROR: Environment '$ENV_NAME' already exists."
    echo "Please choose a different environment name or remove the existing environment manually."
    echo ""
    echo "To remove it manually:"
    echo "    pyenv uninstall '$ENV_NAME'"
    echo ""
    exit 1
fi

# install Python if necessary

if ! pyenv versions --bare | grep -q "^${PYTHON_VERSION}$"; then

    echo ""
    echo "Installing Python $PYTHON_VERSION..."

    pyenv install "$PYTHON_VERSION"
fi

# create environment

echo ""
echo "Creating environment '$ENV_NAME'..."

pyenv virtualenv "$PYTHON_VERSION" "$ENV_NAME"

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create environment."
    exit 1
fi

# get path of new environment

ENV_PATH=$(pyenv prefix "$ENV_NAME")

echo ""
echo "New environment:"
echo "$ENV_PATH"

# use Python and pip directly from the new environment
ENV_PYTHON="$ENV_PATH/bin/python"
ENV_PIP="$ENV_PATH/bin/pip"

echo ""
echo "Python executable:"
echo "$ENV_PYTHON"

echo ""
echo "Pip executable:"
echo "$ENV_PIP"

# extract pip packages from txt

TEMP_REQUIREMENTS=$(mktemp)

awk '
    /^pip:/ {
        in_pip=1
        next
    }

    in_pip && /^  - / {
        sub(/^  - /, "")
        print
        next
    }

    in_pip && !/^  - / {
        exit
    }
' "$ENV_FILE" > "$TEMP_REQUIREMENTS"

echo ""
echo "Packages to install:"
echo "----------------------------------------"
cat "$TEMP_REQUIREMENTS"
echo "----------------------------------------"

# install packages into new environment

echo ""
echo "Installing packages into $ENV_NAME..."

"$ENV_PYTHON" -m pip install --upgrade pip

"$ENV_PYTHON" -m pip install -r "$TEMP_REQUIREMENTS"

if [ $? -ne 0 ]; then
    echo ""
    echo "ERROR: Package installation failed."
    rm "$TEMP_REQUIREMENTS"
    exit 1
fi

rm "$TEMP_REQUIREMENTS"

# verify installation

echo ""
echo "========================================"
echo "Environment successfully created!"
echo "========================================"

echo ""
echo "Environment:"
echo "    $ENV_NAME"

echo ""
echo "Python:"
"$ENV_PYTHON" --version

echo ""
echo "Installed packages:"
"$ENV_PYTHON" -m pip freeze

echo ""
echo "To activate it:"
echo "    pyenv activate $ENV_NAME"