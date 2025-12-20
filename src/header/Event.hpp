#ifndef EVENT_H
#define EVENT_H

enum EventType { ARRIVAL, DEPARTURE };

struct Event {
    double time;
    EventType type;

    bool operator>(const Event& other) const {
        return time > other.time;
    }
};

#endif
