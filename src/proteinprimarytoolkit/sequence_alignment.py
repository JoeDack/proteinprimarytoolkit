from time import sleep

from utils import *

BASE_CLUSTAL_API_URL = r"https://www.ebi.ac.uk/Tools/services/rest/clustalo"

def pairwise_sequence_alignment(sequence1: str, 
                                sequence2: str, 
                                match_score: int = 1, 
                                mismatch_score: int = -1, 
                                indel_score: int = -1) -> tuple[str, str]:
    """
    Pairwise sequence alignment using the Needleman-Wunsch algorithm
    """

    length1 = len(sequence1)
    length2 = len(sequence2)

    grid = [[0] * (length1 + 1) for _ in range(length2 + 1)]
    grid[0] = [i * indel_score for i in range(length1 + 1)]
    
    for i in range(length2 + 1):
        grid[i][0] = i * indel_score

    for i in range(1, length2 + 1):
        row = grid[i]
        previous_row = grid[i-1]
        c2 = sequence2[i-1]

        for j in range(1, length1 + 1):
            # Top and left (vertical/horizontal) scores represent an indel
            # Top left (diagonal) represents a match or mismatch
            # Choose the highest score out of each one
            left = row[j-1]
            top = previous_row[j]
            top_left = previous_row[j-1]
            c1 = sequence1[j-1]

            diagonal_score = top_left + (match_score if c1 == c2 else mismatch_score)
            grid[i][j] = max(diagonal_score, top + indel_score, left + indel_score)
        
    
    aligned1 = []
    aligned2 = []

    i = length2
    j = length1

    # Tracing back to the top left
    while i > 0 or j > 0:
        # Prefer diagonal when possible
        if i > 0 and j > 0:
            if sequence1[j-1] == sequence2[i-1]:
                score = match_score
            else:
                score = mismatch_score

            if grid[i][j] == grid[i-1][j-1] + score:
                aligned1.append(sequence1[j-1])
                aligned2.append(sequence2[i-1])

                i -= 1
                j -= 1

                continue
        
        if i > 0 and grid[i][j] == grid[i-1][j] + indel_score:
            aligned1.append("-")
            aligned2.append(sequence2[i-1])
            
            i -= 1

            continue

        if j > 0 and grid[i][j] == grid[i][j-1] + indel_score:
            aligned1.append(sequence1[j-1])
            aligned2.append("-")

            j -= 1

            continue
    
    aligned1.reverse()
    aligned2.reverse()

    return "".join(aligned1), "".join(aligned2)


def multi_sequence_alignment(sequences: str, sequence_type: str, email: str, timeout: int | float = DEFAULT_TIMEOUT, wait_length: int | float = 0.1) -> str:
    """
    Multi-sequence alignment using the Clustal Omega API. Sequences should be in FASTA format.
    """
    session = get_session()

    data = {"email": email,
            "sequence": sequences,
            "stype": sequence_type
            }

    response = session.post(f"{BASE_CLUSTAL_API_URL}/run", data=data, timeout=timeout)
    raise_for_api_status(response)
    # API returns only job ID
    job_id = response.text.strip()


    status_url = fr"{BASE_CLUSTAL_API_URL}/status/{job_id}"
    status = ""

    # Job takes time - wait until job finishes before requesting result
    while True:
        status = get_text_from_api(session, status_url, timeout).strip()

        if status == "FINISHED":
            break
        elif status in {"QUEUED", "RUNNING"}:
            # Wait a bit to avoid bombarding the server with requests
            sleep(wait_length)
        elif status in {"ERROR", "FAILURE"}:
            raise APIError(f"Clustal Omega job failed with status: {status}")
        elif status == "NOT_FOUND":
            raise APIError(f"Clustal Omega job not found or expired")
        else:
            raise APIError(f"Unexpected status code: {status}")

    
    result_url = fr"{BASE_CLUSTAL_API_URL}/result/{job_id}/aln-clustal"
    result = get_text_from_api(session, result_url, timeout)

    return result