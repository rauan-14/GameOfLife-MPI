#include <iostream>
#include <mpi.h>
#include <vector>
#include <algorithm>
using namespace std;

struct PData {
	int p_id;
	int p_row, p_col;
	int num_rows, num_cols;
	char** board;
	bool border_state[8] = {false, false, false, false, false, false, false, false};
	char* send_buffer;
	char* recv_buffer;
	/* Index mapping */
	/*
	0 - left
	1 - top-left
	2 - top
	3 - top-right
	4 - right
	5 - bottom-right
	6 - bottom
	7 - bottom-left
	*/
	/* Cell status */
	/*
	'0' - dead
	'1' - alive
	*/

	/* Update util functions */
	void updateCell(char data, int i, int j) {
		if(data > '0') {
			if(board[i][j] > '0') {
				board[i][j]++;
			}
			else {
				board[i][j]--;
			}
		}
	}
	void updateRow(char* data, int row) {
		for(int i = 0; i < num_cols; i++) {
			for(int j = i-1; j <= i+1; j++) {
				if(j >= 0 && j < num_cols)
					updateCell(data[i], row, j);
			}
		}
	}
	void updateCol(char* data, int col) {
		for(int i = 0; i < num_rows; i++) {
			for(int j = i-1; j <= i+1; j++) {
				if(j >= 0 && j < num_rows)
					updateCell(data[i], j, col);
			}
		}
	}

	void updateInner() {
		for(int i = 0; i < num_rows; i++) {
			for(int j = 0; j < num_cols; j++) {
				for(int r = i-1; r <= i+1; r++) {
					for(int c = j-1; c <= j+1; c++) {
						if(i == r && j == c)
							continue;
						if(r >= 0 && r < num_rows && c >= 0 && c < num_cols) {
							updateCell(board[i][j], r, c);
						}
					}
				}
			}
		}
	}
	void beautify() {
		for(int i = 0; i < num_rows; i++) {
			for(int j = 0; j < num_cols; j++) {
				if(board[i][j]-'0' == -3 || board[i][j]-'0' == 3 || board[i][j]-'0' == 4)
					board[i][j] = '1';
				else
					board[i][j] = '0';
			}
		}
	}
	void updateLeft(char *data) {
		updateCol(data, 0);
	}
	void updateTopLeft(char* data) {
		updateCell(data[0], 0, 0);
	}
	void updateTop(char* data) {
		updateRow(data, 0);
	}
	void updateTopRight(char* data) {
		updateCell(data[0], 0, num_cols-1);
	}
	void updateRight(char* data) {
		updateCol(data, num_cols-1);
	}
	void updateBottomRight(char* data) {
		updateCell(data[0], num_rows-1, num_cols-1);
	}
	void updateBottom(char* data) {
		updateRow(data, num_rows-1);
	}
	void updateBottomLeft(char* data) {
		updateCell(data[0], num_rows-1, 0);
	}
};

int processNum(int i, int j, int size, int dims[]) {
	int row_chunk = size / dims[0];
	int col_chunk = size / dims[1];
	int p_row = i / row_chunk;
	int p_col = j / col_chunk;
	if(p_row >= dims[0])
		p_row = dims[0] - 1;
	if(p_col >= dims[1])
		p_col = dims[1] - 1;
	int rank = dims[1] * p_row + p_col;
	return rank;
}

pair<int, int> getProcessCoordinates(int p_id, int size, int dims[]) {
	int p_row = p_id / dims[1];
	int p_col = p_id % dims[1];
	return make_pair(p_row, p_col);
}

pair<int, int> getProcessBoardSize(int p_id, int size, int dims[]) {

	int num_rows = size / dims[0];
	pair<int, int> coord = getProcessCoordinates(p_id, size, dims);
	if(coord.first == dims[0] - 1 && size % dims[0] != 0)
		num_rows += size % dims[0];
	
	int num_cols = size / dims[1];
	if(coord.second == dims[1] - 1 && size % dims[1] != 0)
		num_cols += size % dims[1];

	return make_pair(num_rows, num_cols);
}


int main(int argc, char** argv) {


	MPI_Init(&argc, &argv);
	
	/* Declare variables */
	int size, num_gens, num_ghost, num_p;
	PData p_data;
	MPI_Comm_rank(MPI_COMM_WORLD, &p_data.p_id);
	MPI_Comm_size(MPI_COMM_WORLD, &num_p);
	
	/* Read in the parameters */
	if(p_data.p_id == 0) {
		int buffer[3];
		cin >> size >> num_gens >> num_ghost;
		buffer[0] = size;
		buffer[1] = num_gens; 
		buffer[2] = num_ghost;
		for(int dest = 1; dest < num_p; dest++) {
			MPI_Send(buffer, 3, MPI_INT, dest, 0, MPI_COMM_WORLD);
		}

	}
	else {
		int buffer[3];
		MPI_Recv(buffer, 3, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		size = buffer[0];
		num_gens = buffer[1];
		num_ghost = buffer[2];
	}
	MPI_Barrier(MPI_COMM_WORLD);

	/* Divide processes into 2D grid */
	int dims[] = {0, 0};
	int periods[] = {0, 0};
	MPI_Dims_create(num_p, 2, dims);

	/* Handle the case when the are more than needed processors */
	dims[0] = min(dims[0], size);
	dims[1] = min(dims[1], size);
	num_p = dims[0] * dims[1];
	/* Create a new communicator */
	MPI_Comm active_comm;
	bool process_active = p_data.p_id < num_p;
	MPI_Comm_split(MPI_COMM_WORLD, process_active, p_data.p_id, &active_comm);
	

	if(p_data.p_id < num_p) {
		/* Calculate needed sizes */
		p_data.p_row = p_data.p_id / dims[1];
		p_data.p_col = p_data.p_id % dims[1];
		
		p_data.num_rows = size / dims[0];
		if(p_data.p_row == dims[0] - 1 && size % dims[0] != 0)
			p_data.num_rows += size % dims[0];
		
		p_data.num_cols = size / dims[1];
		if(p_data.p_col == dims[1] - 1 && size % dims[1] != 0)
			p_data.num_cols += size % dims[1];


		/* Dynamically allocate the board */
		p_data.board = new char*[p_data.num_rows];
		for(int i = 0; i < p_data.num_rows; i++)
			p_data.board[i] = new char[p_data.num_cols];
		MPI_Barrier(active_comm);
		
		/* Read in the data */
		int max_size = (size / dims[0] + dims[0] - 1) * (size / dims[1] + dims[1] -1);
		if(p_data.p_id == 0) {
			int* indices = new int[num_p];
			char** data = new char*[num_p];
			for(int i = 0; i < num_p; i++) {
				indices[i] = 0;
				data[i] = new char[max_size];
			}
			char cell;
			for(int i = 0; i < size; i++) {
				for(int j = 0; j < size; j++) {
					int id = processNum(i , j, size, dims);
					cin >> cell;
					cell = cell == '#' ? '1' : '0';
					data[id][indices[id]++] = cell;
				}
			}
			for(int i = 1; i < num_p; i++) {
				MPI_Send(data[i], indices[i], MPI_CHAR, i, 0, MPI_COMM_WORLD);
			}
			for(int i = 0; i < p_data.num_rows; i++) {
				for(int j = 0; j < p_data.num_cols; j++) {
					p_data.board[i][j] = data[0][i * p_data.num_cols + j];
				}
			}
			for(int i = 0; i < num_p; i++)
				delete[] data[i];
			delete[] data;
		}
		else {
			char* data = new char[max_size];
			MPI_Recv(data, p_data.num_rows * p_data.num_cols, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for(int i = 0; i < p_data.num_rows; i++) {
				for(int j = 0; j < p_data.num_cols; j++) {
					p_data.board[i][j] = data[i * p_data.num_cols + j];
				}
			}
			delete[] data;
		}


		/* Update border states */
		/* Left */
		if(p_data.p_col > 0) {
			p_data.border_state[0] = true;
		}
		/* Top-left */
		if(p_data.p_col > 0 && p_data.p_row > 0) {
			p_data.border_state[1] = true;
		}
		/* Top */
		if(p_data.p_row > 0) {
			p_data.border_state[2] = true;
		}
		/* Top-right */
		if(p_data.p_row > 0 && p_data.p_col < dims[1]-1) {
			p_data.border_state[3] = true;
		}
		/* Right */
		if(p_data.p_col < dims[1]-1) {
			p_data.border_state[4] = true;
		}
		/* Bottom-right */
		if(p_data.p_row < dims[0]-1 && p_data.p_col < dims[1]-1) {
			p_data.border_state[5] = true;
		}
		/* Bottom */
		if(p_data.p_row < dims[0]-1) {
			p_data.border_state[6] = true;
		}
		/* Bottom-left */
		if(p_data.p_row < dims[0]-1 && p_data.p_col > 0) {
			p_data.border_state[7] = true;
		}

		/* Prepare send and receive buffers */
		max_size = max(p_data.num_cols, p_data.num_rows);
		p_data.send_buffer = new char[max_size];
		p_data.recv_buffer = new char[max_size];
		MPI_Request request;

		/* Start simulating */
		for(int gen = 0; gen < num_gens; gen++) {
			/* Send left */
			if(p_data.border_state[0]) {
				for(int i = 0; i < p_data.num_rows; i++)
					p_data.send_buffer[i] = p_data.board[i][0];
				MPI_Isend(p_data.send_buffer, p_data.num_rows, MPI_CHAR, p_data.p_id-1, 0, MPI_COMM_WORLD, &request);
			}
			if(p_data.border_state[4]) {
				MPI_Recv(p_data.recv_buffer, p_data.num_rows, MPI_CHAR, p_data.p_id+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				p_data.updateRight(p_data.recv_buffer);
			}
			MPI_Barrier(active_comm);

			/* Send top-left */
			if(p_data.border_state[1]) {
				p_data.send_buffer[0] = p_data.board[0][0];
				MPI_Isend(p_data.send_buffer, 1, MPI_CHAR, p_data.p_id-dims[1]-1, 0, MPI_COMM_WORLD, &request);
			}
			if(p_data.border_state[5]) {
				MPI_Recv(p_data.recv_buffer, 1, MPI_CHAR, p_data.p_id+dims[1]+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				p_data.updateBottomRight(p_data.recv_buffer);
			}
			MPI_Barrier(active_comm);

			/* Send top */
			if(p_data.border_state[2]) {
				MPI_Isend(p_data.board[0], p_data.num_cols, MPI_CHAR, p_data.p_id-dims[1], 0, MPI_COMM_WORLD, &request);
			}
			if(p_data.border_state[6]) {
				MPI_Recv(p_data.recv_buffer, p_data.num_cols, MPI_CHAR, p_data.p_id+dims[1], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				p_data.updateBottom(p_data.recv_buffer);
			}
			MPI_Barrier(active_comm);

			/* Send top-right */
			if(p_data.border_state[3]) {
				p_data.send_buffer[0] = p_data.board[0][p_data.num_cols-1];
				MPI_Isend(p_data.send_buffer, 1, MPI_CHAR, p_data.p_id-dims[1]+1, 0, MPI_COMM_WORLD, &request);
			}
			if(p_data.border_state[7]) {
				MPI_Recv(p_data.recv_buffer, 1, MPI_CHAR, p_data.p_id+dims[1]-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				p_data.updateBottomLeft(p_data.recv_buffer);
			}
			MPI_Barrier(active_comm);

			/* Send right */
			if(p_data.border_state[4]) {
				for(int i = 0; i < p_data.num_rows; i++)
					p_data.send_buffer[i] = p_data.board[i][p_data.num_cols-1];
				MPI_Isend(p_data.send_buffer, p_data.num_rows, MPI_CHAR, p_data.p_id+1, 0, MPI_COMM_WORLD, &request);
			}
			if(p_data.border_state[0]) {
				MPI_Recv(p_data.recv_buffer, p_data.num_rows, MPI_CHAR, p_data.p_id-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				p_data.updateLeft(p_data.recv_buffer);
			}
			MPI_Barrier(active_comm);

			/* Send bottom-right */
			if(p_data.border_state[5]) {
				p_data.send_buffer[0] = p_data.board[p_data.num_rows-1][p_data.num_cols-1];
				MPI_Isend(p_data.send_buffer, 1, MPI_CHAR, p_data.p_id+dims[1]+1, 0, MPI_COMM_WORLD, &request);
			}
			if(p_data.border_state[1]) {
				MPI_Recv(p_data.recv_buffer, 1, MPI_CHAR, p_data.p_id-dims[1]-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				p_data.updateTopLeft(p_data.recv_buffer);
			}

			MPI_Barrier(active_comm);

			/* Send bottom */
			if(p_data.border_state[6]) {
				MPI_Isend(p_data.board[p_data.num_rows-1], p_data.num_cols, MPI_CHAR, p_data.p_id+dims[1], 0, MPI_COMM_WORLD, &request);
			}
			if(p_data.border_state[2]) {
				MPI_Recv(p_data.recv_buffer, p_data.num_cols, MPI_CHAR, p_data.p_id-dims[1], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				p_data.updateTop(p_data.recv_buffer);
			}

			MPI_Barrier(active_comm);

			/* Send bottom-left */
			if(p_data.border_state[7]) {
				p_data.send_buffer[0] = p_data.board[p_data.num_rows-1][0];
				MPI_Isend(p_data.send_buffer, 1, MPI_CHAR, p_data.p_id+dims[1]-1, 0, MPI_COMM_WORLD, &request);
			}
			if(p_data.border_state[3]) {
				MPI_Recv(p_data.recv_buffer, 1, MPI_CHAR, p_data.p_id-dims[1]+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				p_data.updateTopRight(p_data.recv_buffer);
			}
			MPI_Barrier(active_comm);

			p_data.updateInner();
			p_data.beautify();
			MPI_Barrier(active_comm);
		}


		/* Gather the data at processor rank 0 and output */
		if(p_data.p_id == 0) {
			char** matrix = new char*[size];
			for(int i = 0; i < size; i++) {
				matrix[i] = new char[size];
			}
			max_size = (size / dims[0] + dims[0] - 1) * (size / dims[1] + dims[1] -1);
			char* recv_buffer = new char[max_size];
			

			int row_chunk = size / dims[0];
			int col_chunk = size / dims[1];

			for(int i = 1; i < num_p; i++) {
				pair<int, int> p_coord = getProcessCoordinates(i, size, dims);
				pair<int, int> p_sizes = getProcessBoardSize(i, size, dims);
				MPI_Recv(recv_buffer, p_sizes.first * p_sizes.second, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				int row_offset = row_chunk * p_coord.first;
				int col_offset = col_chunk * p_coord.second;
				for(int i = 0; i < p_sizes.first; i++) {
					for(int j = 0; j < p_sizes.second; j++) {
						matrix[row_offset+i][col_offset+j] = recv_buffer[i*p_sizes.second+j];
					}
				}
			}
			for(int i = 0; i < p_data.num_rows; i++) {
				for(int j = 0; j < p_data.num_cols; j++) {
					matrix[i][j] = p_data.board[i][j];
				}
			}
			for(int i = 0; i < size; i++) {
				for(int j = 0; j < size; j++) {
					if(matrix[i][j] == '1')
						cout << "#";
					else
						cout << '.';
				}
				cout << endl;
			}

		}

		else {
			char* data = new char[p_data.num_rows * p_data.num_cols];
			for(int i = 0; i < p_data.num_rows; i++) {
				for(int j = 0; j < p_data.num_cols; j++) {
					data[i * p_data.num_cols + j] = p_data.board[i][j];
				}
			}
			MPI_Send(data, p_data.num_rows * p_data.num_cols, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
}

