function format-and-lint {  [[ -z $1 ]] && { echo "Provide target directory" >&2; exit 1; }

  source venv/bin/activate
  python -m isort --profile black $1
  python -m black --line-length 120 $1
  python -m flake8 --ignore=E203,W503 --max-line-length 120 --statistics $1
}