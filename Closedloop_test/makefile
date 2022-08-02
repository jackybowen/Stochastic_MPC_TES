IMG_NAME=julia

COMMAND_RUN=docker run \
	  --name ${IMG_NAME} \
	  --detach=false \
	  --rm \
	  -i \
	  -t \
	  -p 127.0.0.1:5000:5000 \
	  ${IMG_NAME} bash

build:
	docker build --network host --no-cache --rm -t ${IMG_NAME} .

remove-image:
	docker rmi ${IMG_NAME}

run:
	$(COMMAND_RUN)
			
