MODULE ?= python

container:
	sudo docker build -f $(MODULE).Dockerfile -t nakor/coconet-paper-$(MODULE) . \
	&& sudo docker push nakor/coconet-paper-$(MODULE) \
	&& rm -f ~/.singularity_images.cache/nakor-coconet-paper-$(MODULE).img

all:
	make container MODULE=python
	make container MODULE=station-aloha
	make container MODULE=camisim
