
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

osascript 2>/dev/null <<EOF
    tell application "System Events"
        tell process "Terminal" to keystroke "t" using command down
    end
    tell application "Terminal"
        activate
        do script with command "cd \"$DIR/\"; $*" in window 0
        do script with command "./osx_check_postgres_redis.sh" in window 0
    end
EOF
sleep 2

osascript 2>/dev/null <<EOF
    tell application "System Events"
        tell process "Terminal" to keystroke "t" using command down
    end
    tell application "Terminal"
        activate
        do script with command "cd \"$DIR/../client\"; $*" in window 0
        do script with command "gulp watch" in window 0
    end
    tell application "System Events"
        tell process "Terminal" to keystroke "t" using command down
    end
    tell application "Terminal"
        activate
        do script with command "cd \"$DIR/..\"; $*" in window 0
        do script with command "python bin/run_server.py" in window 0
    end
    tell application "System Events"
        tell process "Terminal" to keystroke "t" using command down
    end
    tell application "Terminal"
        activate
        do script with command "cd \"$DIR/..\"; $*" in window 0
        do script with command "celery -A server.webapp.tasks.celery_instance worker -l info" in window 0
    end
EOF
