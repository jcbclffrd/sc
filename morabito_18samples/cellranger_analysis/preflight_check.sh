#!/bin/bash

# Pre-flight checklist for WSL CellRanger setup
# Run this on the ARM machine BEFORE starting WSL setup

echo "============================================================"
echo "PRE-FLIGHT CHECKLIST FOR WSL CELLRANGER SETUP"
echo "============================================================"
echo ""

# Color codes
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

checks_passed=0
checks_failed=0

# Check 1: FASTQ files exist
echo "1. Checking FASTQ files..."
if [ -d "/home/jacobc/sc/sra_downloads" ]; then
    fastq_count=$(ls /home/jacobc/sc/sra_downloads/*.fastq.gz 2>/dev/null | wc -l)
    if [ "$fastq_count" -eq 36 ]; then
        echo -e "   ${GREEN}✓ Found 36 FASTQ files (18 samples × 2 files)${NC}"
        ((checks_passed++))
    else
        echo -e "   ${RED}✗ Found $fastq_count FASTQ files (expected 36)${NC}"
        ((checks_failed++))
    fi
else
    echo -e "   ${RED}✗ Directory /home/jacobc/sc/sra_downloads not found${NC}"
    ((checks_failed++))
fi

# Check 2: Sample mapping exists
echo "2. Checking sample mapping file..."
if [ -f "/home/jacobc/sc/morabito_18samples/sample_mapping.csv" ]; then
    sample_count=$(tail -n +2 /home/jacobc/sc/morabito_18samples/sample_mapping.csv | wc -l)
    echo -e "   ${GREEN}✓ Found sample_mapping.csv with $sample_count samples${NC}"
    ((checks_passed++))
else
    echo -e "   ${RED}✗ sample_mapping.csv not found${NC}"
    ((checks_failed++))
fi

# Check 3: STAR merged file exists (for comparison)
echo "3. Checking STAR merged output..."
if [ -f "/home/jacobc/sc/morabito_18samples/merged_18samples_genes.h5ad" ]; then
    size=$(ls -lh /home/jacobc/sc/morabito_18samples/merged_18samples_genes.h5ad | awk '{print $5}')
    echo -e "   ${GREEN}✓ Found merged_18samples_genes.h5ad ($size)${NC}"
    ((checks_passed++))
else
    echo -e "   ${RED}✗ merged_18samples_genes.h5ad not found${NC}"
    ((checks_failed++))
fi

# Check 4: SSH server running
echo "4. Checking SSH server..."
if systemctl is-active --quiet ssh || systemctl is-active --quiet sshd; then
    echo -e "   ${GREEN}✓ SSH server is running${NC}"
    ((checks_passed++))
else
    echo -e "   ${YELLOW}⚠ SSH server may not be running${NC}"
    echo "   To enable: sudo systemctl start ssh"
    ((checks_failed++))
fi

# Check 5: Get IP address
echo "5. Network configuration..."
ip_addr=$(hostname -I | awk '{print $1}')
if [ -n "$ip_addr" ]; then
    echo -e "   ${GREEN}✓ IP Address: $ip_addr${NC}"
    echo "   Use this IP for WSL sshfs mount"
    ((checks_passed++))
else
    echo -e "   ${RED}✗ Could not determine IP address${NC}"
    ((checks_failed++))
fi

# Check 6: Repository status
echo "6. Checking git repository..."
cd /home/jacobc/sc
if [ -d ".git" ]; then
    branch=$(git branch --show-current)
    status=$(git status --porcelain | wc -l)
    echo -e "   ${GREEN}✓ Git repository initialized (branch: $branch)${NC}"
    if [ "$status" -gt 0 ]; then
        echo -e "   ${YELLOW}⚠ $status uncommitted changes${NC}"
        echo "   Consider: git add cellranger_analysis/ && git commit -m 'Add CellRanger setup'"
    fi
    ((checks_passed++))
else
    echo -e "   ${YELLOW}⚠ Not a git repository${NC}"
    echo "   Consider: git init && git add . && git commit -m 'Initial commit'"
fi

# Check 7: CellRanger analysis directory
echo "7. Checking CellRanger analysis directory..."
if [ -d "/home/jacobc/sc/morabito_18samples/cellranger_analysis" ]; then
    file_count=$(ls /home/jacobc/sc/morabito_18samples/cellranger_analysis/*.{sh,py,md} 2>/dev/null | wc -l)
    echo -e "   ${GREEN}✓ Found cellranger_analysis directory with $file_count files${NC}"
    ((checks_passed++))
else
    echo -e "   ${RED}✗ cellranger_analysis directory not found${NC}"
    ((checks_failed++))
fi

# Summary
echo ""
echo "============================================================"
echo "SUMMARY"
echo "============================================================"
echo -e "Checks passed: ${GREEN}$checks_passed${NC}"
echo -e "Checks failed: ${RED}$checks_failed${NC}"
echo ""

if [ "$checks_failed" -eq 0 ]; then
    echo -e "${GREEN}✓ All checks passed! Ready for WSL setup.${NC}"
    echo ""
    echo "Next steps:"
    echo "1. Note the IP address: $ip_addr"
    echo "2. Verify SSH access from WSL: ssh jacobc@$ip_addr"
    echo "3. On WSL machine, follow AGENT_PROMPT.md"
else
    echo -e "${RED}⚠ Some checks failed. Please resolve issues before proceeding.${NC}"
fi

echo ""
echo "============================================================"
echo "INFORMATION FOR WSL SETUP"
echo "============================================================"
echo ""
echo "Copy these details to use on WSL machine:"
echo ""
echo "ARM_MACHINE_IP=$ip_addr"
echo "ARM_USERNAME=jacobc"
echo "FASTQ_PATH=/home/jacobc/sc/sra_downloads"
echo "REPO_PATH=/home/jacobc/sc"
echo ""
echo "Test SSH connection from WSL:"
echo "  ssh jacobc@$ip_addr 'ls /home/jacobc/sc/sra_downloads/*.fastq.gz | wc -l'"
echo "  (Should return: 36)"
echo ""
echo "============================================================"
