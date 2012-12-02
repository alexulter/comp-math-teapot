cd png
mencoder "mf://*.png" -mf fps=5 -o ../output.avi -ovc lavc -lavcopts vcodec=mpeg4:vbitrate=6000
cd ..
