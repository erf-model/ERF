#ifndef ERF_PLOTFILE_STREAM_H_
#define ERF_PLOTFILE_STREAM_H_

namespace erf_plotfile
{

template <typename T>
constexpr T stream_format (const int stream, const T stream_1, const T stream_2) noexcept
{
    return stream == 2 ? stream_2 : stream_1;
}

} // namespace erf_plotfile

#endif
