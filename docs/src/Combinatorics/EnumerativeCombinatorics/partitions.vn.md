```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Các phân hoạch (partitions)

Một **phân hoạch** (partition) của một số nguyên không âm $n$ là một dãy giảm dần (decreasing sequence) $\lambda_1 \geq \lambda_2 \geq \dots \geq \lambda_r$ của các số nguyên dương $\lambda_i$ sao cho $n = \lambda_1 + \dots + \lambda_r$.
Các $\lambda_i$ được gọi là các **phần tử** (part) của phân hoạch và $r$ được gọi là **độ dài** (length).
Tham khảo thêm tại [Ful97](@cite) và [Knu11](@cite), Mục 7.2.1.4.

Một phân hoạch có thể được mã hóa (encoded) dưới dạng một mảng (array) với các phần tử (element) là $\lambda_i$.
OSCAR cung cấp kiểu tham số (parametric type) `Partition{T}` và đây là một kiểu con (subtype) của `AbstractVector{T}`.
Ở đây, `T` có thể là bất kỳ kiểu con nào của `IntegerUnion`.
Hiệu suất không bị hảnh hưởng khi sử dụng kiểu riêng cho phân hoạch thay vì chỉ sử dụng mảng.
Kiểu tham số này cho phép tăng hiệu suất bằng cách sử dụng các kiểu số nguyên nhỏ hơn.

```@docs
partition
```
Bởi vì `Partition` là một kiểu con của `AbstractVector`, nên tất cả các hàm có thể được sử dụng cho vector (mảng một chiều) cũng có thể được sử dụng cho các phân hoạch.
```jldoctest
julia> P = partition(6, 4, 4, 2)
[6, 4, 4, 2]

julia> length(P)
4

julia> P[1]
6
```
Tuy nhiên, thông thường, $|\lambda| := n$ được gọi là **kích cỡ** (size) của $\lambda$.
Trong Julia, hàm `size` cho mảng đã tồn tại và trả về *kích cỡ* của mảng.
Thay vào đó, bạn có thể sử dụng hàm `sum` của Julia để lấy tổng của các phần tử.
```jldoctest
julia> P = partition(6, 4, 4, 2)
[6, 4, 4, 2]

julia> sum(P)
16
```

Trong các thuật toán liên quan đến phân hoạch, sẽ có lúc ta truy cập các phần tử vượt quá độ dài của phân hoạch. Để thuận tiện, ta thường kì vọng nhận giá trị bằng 0 thay vì bị bắt lỗi.
Vì vậy, OSCAR cung cấp hàm `getindex_safe`:
```@docs
getindex_safe
```
Nếu bạn chắc chắn rằng `P[i]` tồn tại, hãy sử dụng `getindex` vì hàm này sẽ nhanh hơn.

## Sinh (generate) và đếm (count) phân hoạch

```@docs
partitions(::Oscar.IntegerUnion)
number_of_partitions(::Oscar.IntegerUnion)
```
Để đếm phân hoạch, ta sử dụng công thức Hardy-Ramanujan-Rademachen, chi tiết [Joh12](@cite).
Xem thêm [Knu11](@cite), Mục 7.2.1.4 và [OEIS](@cite), [A000041](https://oeis.org/A000041).

### Phân hoạch với các hạn chế (restrictions)

```@docs
> Có bao nhiêu cách để trả một euro với điều kiện chỉ sử dụng các đồng xu có giá trị 1, 2, 5, 10, 20, 50, và/hoặc 100 xu?
> Điều gì xảy ra nếu bạn chỉ được phép sử dụng tối đa hai đồng xu của mỗi loại?

Đây là Bài tập 11 trong [Knu11](@cite), Mục 7.2.1.4. Nó bắt nguồn từ bài toán nổi tiếng
"Cách đổi một đô la", xem [Pol56](@cite). Nói chung, vấn đề là tạo ra và/hoặc đếm các phân hoạch thỏa mãn một số hạn chế nhất định. Tất nhiên, bạn có thể tạo ra danh sách tất cả các phân hoạch của 100 (có khoảng 190 triệu) và sau đó lọc kết quả
theo các hạn chế. Nhưng đối với một số loại hạn chế nhất định, ta sẽ có các thuật toán hiệu quả hơn nhiều.
Các hàm trong phần này thực hiện một số trong số đó. Kết hợp với hàm [filter](https://docs.julialang.org/en/v1/base/collections/#Base.filter) của Julia,
bạn cũng có thể xử lý các loại hạn chế tổng quát hơn.

Ví dụ, có chính xác 6 cách cho câu hỏi thứ hai trong bài tập được trích dẫn ở trên:
```jldoctest
julia> collect(partitions(100, [1, 2, 5, 10, 20, 50], [2, 2, 2, 2, 2, 2]))
6-element Vector{Partition{Int64}}:
 [50, 50]
 [50, 20, 20, 10]
 [50, 20, 20, 5, 5]
 [50, 20, 10, 10, 5, 5]
 [50, 20, 20, 5, 2, 2, 1]
 [50, 20, 10, 10, 5, 2, 2, 1]
```
và có 4562 cách cho câu hỏi thứ nhất trong bài tập:
```jldoctest
julia> length(collect(partitions(100, [1, 2, 5, 10, 20, 50])))
4562
```
Bài toán gốc "Cách đổi một đô la" có 292 cách giải:
```jldoctest
julia> length(collect(partitions(100, [1, 5, 10, 25, 50])))
292
```
```@docs
number_of_partitions(::Oscar.IntegerUnion, ::Oscar.IntegerUnion)
```
Để đếm phân hoạch, ta sử dụng quan hệ hồi quy (recurrence relation) $p_k(n) = p_{k - 1}(n - 1) + p_k(n - k)$, trong đó $p_k(n)$ biểu thị số lượng phân hoạch của $n$ thành $k$ phần tử (part); xem [Knu11](@cite), Mục 7.2.1.4, Phương trình (39), và cũng xem [OEIS](@cite), [A008284](https://oeis.org/A008284).

```@docs
partitions(::Oscar.IntegerUnion, ::Oscar.IntegerUnion, ::Oscar.IntegerUnion, ::Oscar.IntegerUnion)
partitions(::T, ::Vector{T}) where T <: Oscar.IntegerUnion
```

## Các phép toán (operations)

*Liên hợp* (conjugate) của một phân hoạch $\lambda$ được tính bằng cách xét biểu đồ Young (Young diagram) của nó
(xem [Tableaux](@ref)) và sau đó lật nó dọc theo đường chéo chính, xem [Ful97](@cite), trang 2, và [Knu11](@cite), Mục 7.2.1.4.
```@docs
conjugate
```

## Các Quan hệ (relations)

Thứ tự **thống trị** (dominance order) trên các phân hoạch là thứ tự từng phần (partial order) $\trianglerighteq$ được định nghĩa bởi $\lambda \trianglerighteq\mu$ nếu và chỉ nếu $\lambda_1 + \dots + \lambda_i \geq \mu_1 + \dots + \mu_i$ cho mọi $i$.
Nếu $\lambda\trianglerighteq\mu$, người ta nói rằng $\lambda$ **thống trị** (dominate) $\mu$.
Xem [Ful97](@cite), trang 26, và [Knu11](@cite), Mục 7.2.1.4, Bài tập 54.

Lưu ý rằng trong khi thứ tự từ điển (lexicographic oder) là một thứ tự toàn phần (total ordering), thứ tự thống trị không phải như vậy.
Hơn nữa, [Knu11](@cite) nói rằng **áp đảo** (majorize) thay vì **thống trị** (dominate) và sử dụng ký hiệu $\succeq$ thay vì $\trianglerighteq$.
```@docs
dominates
```
