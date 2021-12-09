aws-genomics-workflows-cn
---

# TODO list

- 存储方面
	- 使用s3fs/goofys/minfs挂载对象存储作为本地磁盘, 用以顺序读写.

- repo中新功能的引入
	- aws-genomics-workflows
    - amazon-ebs-autoscale

- 标签方面
    - ~~ebs-autoscale的资源贴标签~~

- AWS中国和AWS全球配置上的差异自动配置
    - IAM: arn:aws: <-> arn:aws-cn
    - 域名: 
    	- s3.region.amazonaws.com <->  s3.region.amazonaws.com.cn
    	- ec2.amazonaws.com <-> ec2.amazonaws.com.cn
