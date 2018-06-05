CC = g++
CFLAGS = -std=c++11 -Wall

closest-pair: closestpair.cpp
	$(CC) $(CFLAGS) -o $@ closestpair.cpp

clean:
	rm -f closest-pair
