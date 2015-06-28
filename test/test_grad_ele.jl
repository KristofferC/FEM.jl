import FEM.get_u_idxs
import FEM.get_grad_idxs

@test get_u_dof_idxs(6, 2, 3, 4) == [1, 2, 7, 8, 13, 14, 19, 20, 21, 22, 23, 24]
@test get_grad_dof_idxs(6, 2, 3, 4) == [3, 4, 5, 6, 9, 10, 11, 12, 15, 16, 17, 18]
@test get_grad_idxs_plane(3, 1) == [1, 2, 5, 6, 9, 10]
@test get_grad_idxs_plane(3, 2) == [3, 4, 7, 8, 11, 12]


@test get_u_dof_idxs(6, 2, 3, 2) == [1, 2, 7, 8, 13, 14, 19, 20, 21, 22, 23, 24]
@test get_grad_dof_idxs(6, 2, 3, 2) == [3, 4, 5, 6, 9, 10, 11, 12, 15, 16, 17, 18]
@test get_grad_idxs_plane(3, 1) == [1, 2, 5, 6, 9, 10]
@test get_grad_idxs_plane(3, 2) == [3, 4, 7, 8, 11, 12]
