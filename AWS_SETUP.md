# AWS ATAC-seq Alignment Setup

Fast ATAC-seq alignment on AWS to avoid local download/compute time.

## Cost Estimate
- Instance: c5.9xlarge at $1.53/hour × 2-3 hours = **$3-5**
- Storage: 1TB EBS at $0.10/GB-month × 3 hours = **$0.50**
- Data transfer OUT: 300GB BAMs × $0.09/GB = **$27**
- **Total: ~$30-35** to complete in 2 hours instead of 20 hours locally

## Step 1: Launch EC2 Instance

### Via AWS Console:
1. Go to https://console.aws.amazon.com/ec2/
2. Click "Launch Instance"
3. Configure:
   - **Name**: `atacseq-alignment-star`
   - **AMI**: Ubuntu Server 22.04 LTS (64-bit ARM or x86)
   - **Instance type**: `c5.9xlarge` (36 cores, 72GB RAM)
   - **Key pair**: Use your existing `awsWebsite.pem` or create new
   - **Storage**: 1TB gp3 SSD
   - **Security group**: Allow SSH (port 22) from your IP
4. Click "Launch Instance"

### Via AWS CLI (Alternative):
```bash
aws ec2 run-instances \
  --image-id ami-0e001c9271cf7f3b9 \
  --instance-type c5.9xlarge \
  --key-name awsWebsite \
  --block-device-mappings '[{"DeviceName":"/dev/sda1","Ebs":{"VolumeSize":1000,"VolumeType":"gp3"}}]' \
  --region us-east-1 \
  --tag-specifications 'ResourceType=instance,Tags=[{Key=Name,Value=atacseq-alignment}]'
```

Get the instance IP:
```bash
aws ec2 describe-instances --filters "Name=tag:Name,Values=atacseq-alignment" \
  --query 'Reservations[0].Instances[0].PublicIpAddress' --output text
```

## Step 2: Connect and Setup

Save the IP address, then connect:
```bash
# From your local machine
export AWS_IP="<YOUR-INSTANCE-IP>"
ssh -i ~/.ssh/awsWebsite.pem ubuntu@$AWS_IP
```

## Step 3: Run Setup Script

On the AWS instance, run the automated setup:
```bash
# Download and run setup script
curl -o setup_aws.sh https://raw.githubusercontent.com/jcbclffrd/sc/master/scripts/setup_aws_atacseq.sh
chmod +x setup_aws.sh
./setup_aws.sh
```

Or manually copy the setup script (see `scripts/setup_aws_atacseq.sh`)

## Step 4: Transfer Genome Index

### Option A: Upload from local (if you have fast upload)
```bash
# From your LOCAL machine
rsync -avz -e "ssh -i ~/.ssh/awsWebsite.pem" \
  ~/sc/star_index/ \
  ubuntu@$AWS_IP:~/sc/star_index/
```

### Option B: Rebuild on AWS (30 minutes, one-time)
```bash
# On AWS instance
cd ~/sc
./scripts/build_star_index_aws.sh
```

## Step 5: Run Alignment

On AWS instance:
```bash
cd ~/sc
nohup bash scripts/align_atacseq_aws.sh > logs_atacseq/aws_alignment.log 2>&1 &

# Monitor progress
tail -f logs_atacseq/aws_alignment.log

# Check completed samples
ls atacseq_aligned/SRR*/Aligned.sortedByCoord.out.bam | wc -l
```

## Step 6: Download Results

From your LOCAL machine:
```bash
# Download all BAM files
rsync -avz -e "ssh -i ~/.ssh/awsWebsite.pem" \
  ubuntu@$AWS_IP:~/sc/atacseq_aligned/ \
  ~/sc/atacseq_aligned/

# Or download just completed samples
rsync -avz -e "ssh -i ~/.ssh/awsWebsite.pem" \
  --include='*/' \
  --include='Aligned.sortedByCoord.out.bam' \
  --include='Log.final.out' \
  --exclude='*' \
  ubuntu@$AWS_IP:~/sc/atacseq_aligned/ \
  ~/sc/atacseq_aligned/
```

## Step 7: Cleanup

**IMPORTANT**: Terminate instance when done to stop charges!

```bash
# From local machine
aws ec2 terminate-instances --instance-ids <INSTANCE-ID>

# Or via console: EC2 → Instances → Right-click → Terminate
```

## Timeline

- Setup: 15 minutes
- Genome index upload/build: 30 minutes (one-time)
- Downloads + alignment: 2 hours (9 parallel jobs)
- Download BAMs: 30 minutes
- **Total: ~3 hours active work**

## Monitoring

Check progress:
```bash
# On AWS instance
watch -n 30 'ls atacseq_aligned/SRR*/Aligned.sortedByCoord.out.bam | wc -l'

# Check resource usage
htop

# Check disk space
df -h

# Check network usage
iftop
```

## Troubleshooting

**Out of disk space:**
```bash
# Clean up temp files
rm -rf atacseq_aligned/SRR*/_STARtmp

# Delete FASTQ files after alignment (if downloaded locally)
rm -rf sra_downloads/ATAC-seq/*.fastq*
```

**Alignment failed:**
```bash
# Check logs for specific sample
cat logs_atacseq/SRR14514129_alignment.log

# Re-run just failed samples
./scripts/align_atacseq_aws.sh
```

## Notes

- SRA data is mirrored in us-east-1, so downloads are fast (~129 MB/s)
- c5.9xlarge has 36 cores, we run 9 samples × 4 threads = 36 cores fully utilized
- Expected completion: 32 samples ÷ 9 parallel = 4 batches × 30 min = 2 hours
- BAM files are ~9GB each, so 32 × 9GB = ~288GB to download
