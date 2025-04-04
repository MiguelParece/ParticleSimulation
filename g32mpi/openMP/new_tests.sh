#!/bin/bash

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Function to compare floating point numbers with tolerance
compare_floats() {
    local val1=$1
    local val2=$2
    local tolerance=0.001
    
    # Handle potential syntax errors from bc
    local diff=$(echo "scale=6; sqrt(($val1 - $val2)^2)" | bc 2>/dev/null)
    if [ $? -ne 0 ]; then
        return 1
    fi
    
    local threshold=$(echo "scale=6; $tolerance" | bc)
    
    if (( $(echo "$diff <= $threshold" | bc -l) )); then
        return 0
    else
        return 1
    fi
}

# Test cases array (first few tests for initial debugging)
declare -a tests=(
    "5893 0.05 3 10 10:0.002 0.035:2" #1
    "8555 0.05 3 10 10:0.016 0.049:1" #2
    "12 100 5 10000 10000:76.732 61.943:2209" #3
    "-11 3500 20 500000 10:1984.878 1625.992:35" #4
)

total_tests=0
passed_tests=0

echo "Starting tests..."
echo "----------------"

for test in "${tests[@]}"; do
    # Split test case into input and expected output
    IFS=":" read -r input expected_coords expected_collisions <<< "$test"
    
    total_tests=$((total_tests + 1))
    echo "Test $total_tests: ./parsim-omp $input"
    
    # Run the program and capture output and time
    start_time=$(date +%s.%N)
    output=$(./parsim-omp $input)
    end_time=$(date +%s.%N)
    execution_time=$(echo "$end_time - $start_time" | bc)
    
    if [ $? -ne 0 ]; then
        echo -e "${RED}Test $total_tests Failed: Program crashed${NC}"
        continue
    fi
    
    # Parse the two lines of output separately
    coords=$(echo "$output" | head -n 1)
    collisions=$(echo "$output" | tail -n 1)
    
    # Split expected and actual coordinates
    read -r expected_x expected_y <<< "$expected_coords"
    read -r actual_x actual_y <<< "$coords"
    
    passed=true
    
    # Compare coordinates
    if ! compare_floats "$actual_x" "$expected_x" || ! compare_floats "$actual_y" "$expected_y"; then
        echo -e "${RED}Coordinates mismatch${NC}"
        echo "Expected: $expected_x $expected_y"
        echo "Got     : $actual_x $actual_y"
        passed=false
    fi
    
    # Compare collisions
    if [ "$collisions" != "$expected_collisions" ]; then
        echo -e "${RED}Collisions mismatch${NC}"
        echo "Expected: $expected_collisions"
        echo "Got     : $collisions"
        passed=false
    fi
    
    if [ "$passed" = true ]; then
        echo -e "${GREEN}Test $total_tests Passed${NC}"
        passed_tests=$((passed_tests + 1))
    else
        echo -e "${RED}Test $total_tests Failed${NC}"
    fi
    
    execution_time_min=$(echo "scale=6; $execution_time / 60" | bc)
    printf "Execution time: %.6f minutes\n" $execution_time_min
    echo "----------------"
done

# Print summary
echo "Test Summary:"
echo "Tests passed: $passed_tests/$total_tests"
if [ $passed_tests -eq $total_tests ]; then
    echo -e "${GREEN}All tests passed!${NC}"
else
    echo -e "${RED}Some tests failed${NC}"
fi
