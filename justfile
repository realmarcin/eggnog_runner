default: lint test

install:
	poetry install

build:
	poetry build

lint:
	poetry run ruff check kg_microbe

test:
	poetry run pytest tests --cov=kg_microbe -s --cov-report term-missing
