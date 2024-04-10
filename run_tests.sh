#!/bin/bash

# Run the Clojure file and redirect output to a text file
lein run -m graph-fri > test_outputs.txt

echo "Test outputs have been saved to test_outputs.txt"
