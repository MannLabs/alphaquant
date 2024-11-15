def sort_by_end(items):
    items.sort(key=lambda x: int(x.split('_')[1]))

def sort_by_start(items):
    items.sort(key=lambda x: int(x.split('_')[0]))

def test_sorting_order_preservation():
    # Initial list
    initial_list = ["3_3", "1_3", "2_2", "4_1", "1_2", "3_2", "2_3", "3_1"]
    expected_result = ["1_2", "1_3", "2_2", "2_3", "3_1", "3_2", "3_3", "4_1"]
    print("Initial list:", initial_list)

    # First sort by the end number
    sort_by_end(initial_list)
    print("After sorting by end number:", initial_list)

    # Then sort by the start number
    sort_by_start(initial_list)
    print("After sorting by start number:", initial_list)

    # Check if the final sort is equal to the expected result
    assert initial_list == expected_result, "Final sort is not equal to the expected result."
    print("Test passed: Order preserved for elements with the same start number.")

test_sorting_order_preservation()