#%%
from llama_index.core.tools.tool_spec.base import BaseToolSpec
from datetime import datetime
from typing import List, Tuple, Dict, Any

def working_hours(day: datetime) -> Tuple[datetime, datetime]:
    return (datetime(day.year, day.month, day.day, 8, 0), datetime(day.year, day.month, day.day, 18, 0))

def create_person_calendar(person: str) -> List[Dict[str, any]]:
    calendar = []
    if person == "Person1":
        calendar.extend([
            {"start": datetime(2025, 5, 1, 9, 0), "end": datetime(2025, 5, 1, 10, 0), 
             "title": "User Requirements Gathering", "attendees": ["Person1", "FakePersonA"]},
            {"start": datetime(2025, 5, 1, 14, 0), "end": datetime(2025, 5, 1, 15, 0), 
             "title": "Design Review", "attendees": ["Person1", "FakePersonB"]},
            {"start": datetime(2025, 5, 2, 11, 0), "end": datetime(2025, 5, 2, 12, 0), 
             "title": "Sprint Planning", "attendees": ["Person1", "Person2"]},
            {"start": datetime(2025, 5, 3, 16, 0), "end": datetime(2025, 5, 3, 17, 0), 
             "title": "Team Sync", "attendees": ["Person1", "Person2", "Person3"]},
            {"start": datetime(2025, 5, 3, 10, 0), "end": datetime(2025, 5, 3, 11, 0), 
             "title": "Project Planning", "attendees": ["Person1", "FakePersonF"]}
        ])
    elif person == "Person2":
        calendar.extend([
            {"start": datetime(2025, 5, 1, 9, 0), "end": datetime(2025, 5, 1, 9, 30), 
             "title": "Scrum Standup", "attendees": ["Person2", "Agile Team"]},
            {"start": datetime(2025, 5, 1, 13, 0), "end": datetime(2025, 5, 1, 14, 0), 
             "title": "Client Demo", "attendees": ["Person2", "FakePersonC"]},
            {"start": datetime(2025, 5, 2, 11, 0), "end": datetime(2025, 5, 2, 12, 0), 
             "title": "Sprint Planning", "attendees": ["Person1", "Person2"]},
            {"start": datetime(2025, 5, 3, 16, 0), "end": datetime(2025, 5, 3, 17, 0), 
             "title": "Team Sync", "attendees": ["Person1", "Person2", "Person3"]},
            {"start": datetime(2025, 5, 3, 14, 0), "end": datetime(2025, 5, 3, 15, 0), 
             "title": "Budget Review", "attendees": ["Person2", "FakePersonG"]}
        ])
    elif person == "Person3":
        calendar.extend([
            {"start": datetime(2025, 5, 1, 10, 0), "end": datetime(2025, 5, 1, 11, 0), 
             "title": "Tech Deep Dive", "attendees": ["Person3", "FakePersonD"]},
            {"start": datetime(2025, 5, 2, 14, 0), "end": datetime(2025, 5, 2, 15, 0), 
             "title": "Architecture Review", "attendees": ["Person3", "FakePersonE"]},
            {"start": datetime(2025, 5, 3, 16, 0), "end": datetime(2025, 5, 3, 17, 0), 
             "title": "Team Sync", "attendees": ["Person1", "Person2", "Person3"]},
            {"start": datetime(2025, 5, 3, 9, 0), "end": datetime(2025, 5, 3, 10, 0), 
             "title": "Code Review", "attendees": ["Person3", "FakePersonH"]}
        ])
    return calendar

people = ["Person1", "Person2", "Person3"]
calendar_schedules = {person: create_person_calendar(person) for person in people}
person1_calendar = create_person_calendar("Person1")
person2_calendar = create_person_calendar("Person2")
person3_calendar = create_person_calendar("Person3")
days = [datetime(2025, 5, 1), datetime(2025, 5, 2), datetime(2025, 5, 3)]
working_days = [working_hours(day) for day in days]

#%%
class CalendarToolSpec(BaseToolSpec):
    """A mock calendar tool spec for demonstration purposes."""

    spec_functions = [
        "find_common_free_slots",
        "add_event_to_calendars"
    ]
    
    def __init__(
        self,
        working_days: List[datetime] = working_days,
        calendar_schedules: Dict[str, List[Dict[str, Any]]] = calendar_schedules,
    ):
        self.working_days = working_days
        self.people = list(calendar_schedules.keys())
        self.calendars = list(calendar_schedules.values())
    
    def find_common_free_slots(self):
        """Tool used to find common free slots in calendars between
        Person1, Person2 and Person3. 
        """
        results = dict()
        for day_start, day_end in self.working_days:
            day_calendars = [self.filter_events_by_day(calendar, day_start, day_end) for calendar in calendars]
            common_free = self.find_common_free_slot(day_calendars, day_start, day_end)
            results[day_start.strftime("%Y-%m-%d")] = common_free
        return results
    
    def add_event_to_calendars(
        self,
        start: datetime, 
        end: datetime, 
        title: str, 
        calendars_mapping: Dict[str, List[Dict[str, any]]]) -> str:
        """
        Adds an event to the calendars of specified people if the time slot is free.
        
        Parameters:
            start (datetime): The start time of the event.
            end (datetime): The end time of the event.
            title (str): Title of the event.
            calendars_mapping (Dict[str, List[Dict[str, any]]]): A mapping from person names to their calendar events.
            
        Returns:
            str: Success message if the event was added, otherwise an error message.
        """
        def is_free(calendar: List[Dict[str, any]], start: datetime, end: datetime) -> bool:
            for event in calendar:
                # Check if the new event overlaps with any existing event.
                if not (end <= event["start"] or start >= event["end"]):
                    return False
            return True

        # Verify availability for each affected person.
        for person in self.people:
            cal = calendars_mapping.get(person)
            if cal is None or not is_free(cal, start, end):
                return f"Time slot not available for {person}."
                
        # Build new event and add to each calendar.
        event = {"start": start, "end": end, "title": title, "attendees": self.people.copy()}
        for person in self.people:
            calendars_mapping[person].append(event)
        return("Event added successfully to calendars: " + ", ".join(self.people))
    
    @staticmethod
    def get_free_intervals(
        busy_intervals: List[Dict[str, datetime]],
        day_start: datetime,
        day_end: datetime
    ):
        """
        Given a list of busy intervals (calendar events), computes free time slots 
        within the working day.
        """
        # Sort events by start time
        busy_intervals = sorted(busy_intervals, key=lambda x: x["start"])
        free_intervals: List[Tuple[datetime, datetime]] = []

        # Check for free time before the first event
        if busy_intervals and busy_intervals[0]["start"] > day_start:
            free_intervals.append((day_start, busy_intervals[0]["start"]))

        # Identify free slots between consecutive meetings
        for i in range(len(busy_intervals) - 1):
            current_end = busy_intervals[i]["end"]
            next_start = busy_intervals[i + 1]["start"]
            if current_end < next_start:
                free_intervals.append((current_end, next_start))

        # Check for free time after the last event
        if busy_intervals and busy_intervals[-1]["end"] < day_end:
            free_intervals.append((busy_intervals[-1]["end"], day_end))

        return free_intervals
    
    @staticmethod
    def intersect_intervals(
        intervals1: List[Tuple[datetime, datetime]], 
        intervals2: List[Tuple[datetime, datetime]]
    ) -> List[Tuple[datetime, datetime]]:
        """
        Computes the intersection of two lists of free intervals.
        """
        i, j = 0, 0
        intersections: List[Tuple[datetime, datetime]] = []

        while i < len(intervals1) and j < len(intervals2):
            start1, end1 = intervals1[i]
            start2, end2 = intervals2[j]
            latest_start = max(start1, start2)
            earliest_end = min(end1, end2)
            if latest_start < earliest_end:
                intersections.append((latest_start, earliest_end))
            if end1 < end2:
                i += 1
            else:
                j += 1
        return intersections
    
    def find_common_free_slot(
        self,
        calendars: List[List[Dict[str, datetime]]], 
        day_start: datetime, 
        day_end: datetime) -> List[Tuple[datetime, datetime]]:
        """
        Finds common available time slots for all participants within a single day.
        """
        free_intervals_list: List[List[Tuple[datetime, datetime]]] = []
        for calendar in calendars:
            free_intervals_list.append(self.get_free_intervals(calendar, day_start, day_end))
        common_free: List[Tuple[datetime, datetime]] = free_intervals_list[0]
        for free_intervals in free_intervals_list[1:]:
            common_free = self.intersect_intervals(common_free, free_intervals)
        return common_free

    @staticmethod
    def filter_events_by_day(calendar: List[Dict[str, any]], 
                            day_start: datetime, 
                            day_end: datetime) -> List[Dict[str, any]]:
        """
        Returns events from the calendar that fall within the given day.
        """
        return [event for event in calendar if event["start"] >= day_start and event["start"] < day_end]

    def find_common_free_slot_multi_day(
        self,
        calendars: List[List[Dict[str, any]]], 
        working_days: List[Tuple[datetime, datetime]],
        persons: List[str] = ["Person1", "Person2", "Person3"]
    ) -> List[datetime]:
        """
        Computes common free slots for multiple days and then returns a sorted list 
        of datetime objects (start and end times of free intervals) combined with 
        the start and end times of meetings where all three participants are involved.
        """
        free_boundaries = set()

        # Process each day separately
        for day_start, day_end in working_days:
            # For each calendar, filter events to just those in the day.
            daily_calendars = [self.filter_events_by_day(calendar, day_start, day_end) for calendar in calendars]
            daily_free = self.find_common_free_slot(daily_calendars, day_start, day_end)
            for start, end in daily_free:
                free_boundaries.add(start)
                free_boundaries.add(end)
        
        # Get full meetings (meetings with all three attendees) across all days.
        full_meetings = self.flag_full_meetings(calendars, persons)
        for meeting in full_meetings:
            free_boundaries.add(meeting["start"])
            free_boundaries.add(meeting["end"])
        
        # Return sorted list of datetime objects.
        return sorted(free_boundaries)
    
    @staticmethod
    def ensure_owner_in_meetings(calendar: List[Dict[str, any]], owner: str) -> None:
        """
        Ensures that each meeting in the given calendar includes the calendar owner's name 
        in the attendees list. If not, adds the owner.
        """
        for event in calendar:
            if "attendees" in event:
                if owner not in event["attendees"]:
                    event["attendees"].append(owner)
            else:
                event["attendees"] = [owner]

    @staticmethod
    def flag_full_meetings(
        calendars: List[List[Dict[str, any]]], 
        required_attendees: List[str]
    ) -> List[Dict[str, any]]:
        """
        Flags meetings that include all the required attendees.
        """
        all_events: List[Dict[str, any]] = []
        for calendar in calendars:
            all_events.extend(calendar)
        
        # Deduplicate events based on (start, end, title)
        unique_events: Dict[Tuple[datetime, datetime, str], Dict[str, any]] = {}
        for event in all_events:
            key = (event["start"], event["end"], event["title"])
            if key not in unique_events:
                unique_events[key] = event.copy()  # copy to avoid modifying original
            else:
                unique_events[key]["attendees"] = list(set(unique_events[key]["attendees"]) | set(event["attendees"]))
        
        flagged: List[Dict[str, any]] = []
        for event in unique_events.values():
            if all(attendee in event["attendees"] for attendee in required_attendees):
                flagged.append(event)
        return flagged

    @staticmethod
    def flag_partial_meetings(
        calendars: List[List[Dict[str, any]]], 
        pair: Tuple[str, str]
    ) -> List[Dict[str, any]]:
        """
        Flags meetings that include both persons specified in the pair.
        """
        all_events: List[Dict[str, any]] = []
        for calendar in calendars:
            all_events.extend(calendar)
        
        unique_events: Dict[Tuple[datetime, datetime, str], Dict[str, any]] = {}
        for event in all_events:
            key = (event["start"], event["end"], event["title"])
            if key not in unique_events:
                unique_events[key] = event.copy()
            else:
                unique_events[key]["attendees"] = list(set(unique_events[key]["attendees"]) | set(event["attendees"]))
        
        flagged: List[Dict[str, any]] = []
        for event in unique_events.values():
            if all(person in event["attendees"] for person in pair):
                flagged.append(event)
        return flagged
# %%
