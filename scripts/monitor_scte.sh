#!/bin/bash
# Monitor scTE progress and memory usage

echo "scTE Quantification Monitor"
echo "============================="
echo ""

# Check if running
if ps aux | grep -q "[s]cTE.*-i.*-o"; then
    echo "✓ scTE is RUNNING"
    echo ""
    
    # Show memory usage
    echo "Memory Usage:"
    free -h | grep -E "Mem|Swap"
    echo ""
    
    # Show CPU usage
    echo "CPU Usage:"
    ps aux | grep "[s]cTE.*-i.*-o" | awk '{printf "  PID: %s, CPU: %s%%, MEM: %s%%, TIME: %s\n", $2, $3, $4, $10}'
    echo ""
else
    echo "✗ scTE is NOT running"
    echo ""
fi

# Count completed samples
COMPLETED=$(find ~/sc/scTE_output -name "*.h5ad" 2>/dev/null | wc -l)
echo "Progress: ${COMPLETED} / 150 samples completed"

# Show current sample from log
echo ""
echo "Current sample:"
tail -5 ~/sc/scte_quantification.log 2>/dev/null | grep -E "Processing SRR|Started:|Completed|INFO.*Detect"

# Show disk space
echo ""
echo "Disk Space:"
df -h /home/jacobc/sc | tail -1
