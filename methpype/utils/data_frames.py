__all__ = ['inner_join_data']


def inner_join_data(left_df, right_df, left_on=None, right_on=None):
    """Helper function that performs an inner join on two data frames
    in a strict, consistent manner.

    Arguments:
        left_df {DataFrame} -- First data frame to merge.
        right_df {DataFrame} -- Second data frame to merge.

    Keyword Arguments:
        left_on {string} -- Column of the first data frame to use as the merge column.
            If not provided, the index will be used. (default: {None})
        right_on {string} -- Column of the first data frame to use as the merge column.
            If not provided, the index will be used. (default: {None})

    Returns:
        [DataFrame] -- A new data frame of the merged values.
    """
    left_index = left_on is None
    right_index = right_on is None

    return left_df.merge(
        right_df,
        how='inner',
        left_index=left_index,
        left_on=left_on,
        right_index=right_index,
        right_on=right_on,
        suffixes=(False, False),
    )
