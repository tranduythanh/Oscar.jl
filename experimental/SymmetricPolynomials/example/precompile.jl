using PackageCompiler
using Dates

sysimage_path = "Oscar.so"

if isfile(sysimage_path)
    last_modified = mtime(sysimage_path)
    println("Existing system image found. Last modified: ", Dates.format(unix2datetime(last_modified), "yyyy-mm-dd HH:MM:SS"))
    print("Do you want to recreate the system image? (y/n): ")
    response = lowercase(strip(readline()))
    if response != "y"
        println("Keeping existing system image. Exiting...")
        exit(0)
    end
end

# Create the system image for Oscar only
create_sysimage(
    [:Oscar], 
    sysimage_path=sysimage_path
)

println("System image created successfully!")
