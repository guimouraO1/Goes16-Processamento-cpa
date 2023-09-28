#!/bin/bash
find /home/guimoura/processamento/output/band* -type f -name "*.gif" -exec rm {} \;
find /home/guimoura/processamento/output/band* -type f -name "*.png" -exec rm {} \;
find /home/guimoura/processamento/output/fdcf* -type f -name "*.*" -exec rm {} \;
find /home/guimoura/processamento/output/rrqpef* -type f -name "*.*" -exec rm {} \;
find /home/guimoura/processamento/output/glm* -type f -name "*.*" -exec rm {} \;
find /home/guimoura/processamento/output/truecolor* -type f -name "*.**" -exec rm {} \;