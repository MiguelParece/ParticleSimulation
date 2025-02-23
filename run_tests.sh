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
    "1 2 3 10 1:1.570 0.056:0"
    "1 1 5 100 1:0.786 0.027:0"
    "-10 3 3 100 10:1.733 1.643:2"
    "-50 10000 200 500000 10:5025.384 5303.928:4"
    "1 5000 100 1000000 4:3936.506 131.472:4"
    "1 5000 100 1000000 100:3899.787 156.291:163"
    "1 5000 20 1000000 10:3918.912 143.364:19"
    "1 1000 3 10000 10000:287.788 261.446:31"
    "3 5000 50 1000000 300:3819.032 25.659:469"
    "3 5000 50 1000000 500:3738.436 58.743:804"
    "-1 1000 30 100000 1000:575.878 370.663:1203"
)

total_tests=0
passed_tests=0

echo "Starting tests..."
echo "----------------"

for test in "${tests[@]}"; do
    # Split test case into input and expected output
    IFS=":" read -r input expected_coords expected_collisions <<< "$test"
    
    total_tests=$((total_tests + 1))
    echo "Test $total_tests: ./parsim $input"
    
    # Run the program and capture output and time
    start_time=$(date +%s.%N)
    output=$(./parsim $input)
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
    
    # Print execution time
    printf "Execution time: %.6f seconds\n" $execution_time
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
