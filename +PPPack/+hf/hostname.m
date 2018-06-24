function outName = hostname
%HOSTNAME (ps-utils): return hostname of machine matlab is running on
%   T = HOSTNAME returns the hostname of the current machine.  If FQDN
%   Only the hostname (without domain) is returned.
%
%$Id: hostname.m 125 2008-03-20 20:19:22Z vincent $

% cache check
persistent cached_name
if ~isempty(cached_name), outName=cached_name; return; end

[retval,tName] = system('hostname');

% some error checks
assert(retval == 0, 'Error running hostname');
assert(~any(tName == '.'), 'Dots found in hostname: is it a fqdn?');

outName = deblank(tName);

%load cache
cached_name=outName;


