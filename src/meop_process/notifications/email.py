from __future__ import annotations

from email.message import EmailMessage
import mimetypes
import os
from pathlib import Path
import smtplib
import subprocess

from ..models import EmailNotificationSettings


def _build_message(
    *,
    subject: str,
    body: str,
    recipients: tuple[str, ...],
    from_address: str,
    attachments: tuple[Path, ...],
) -> EmailMessage:
    message = EmailMessage()
    message["Subject"] = subject
    message["To"] = ", ".join(recipients)
    message["From"] = from_address
    message.set_content(body)

    for path in attachments:
        content = path.read_bytes()
        mime_type, _ = mimetypes.guess_type(str(path))
        if mime_type:
            maintype, subtype = mime_type.split("/", 1)
        else:
            maintype, subtype = "application", "octet-stream"
        message.add_attachment(content, maintype=maintype, subtype=subtype, filename=path.name)
    return message


def send_email_message(
    settings: EmailNotificationSettings,
    *,
    subject: str,
    body: str,
    attachments: tuple[Path, ...] = (),
) -> None:
    recipients = tuple(item for item in settings.to if item)
    if not settings.enabled or not recipients:
        return

    from_address = settings.transport.from_address or recipients[0]
    message = _build_message(
        subject=subject,
        body=body,
        recipients=recipients,
        from_address=from_address,
        attachments=attachments,
    )

    transport = settings.transport.transport.strip().lower()
    if transport == "sendmail":
        subprocess.run(
            [settings.transport.sendmail_path, "-t", "-oi"],
            input=message.as_bytes(),
            check=True,
        )
        return

    username = os.getenv(settings.transport.username_env, "") if settings.transport.username_env else ""
    password = os.getenv(settings.transport.password_env, "") if settings.transport.password_env else ""
    if settings.transport.use_ssl:
        with smtplib.SMTP_SSL(settings.transport.host, settings.transport.port) as client:
            if username or password:
                client.login(username, password)
            client.send_message(message)
        return

    with smtplib.SMTP(settings.transport.host, settings.transport.port) as client:
        if settings.transport.starttls:
            client.starttls()
        if username or password:
            client.login(username, password)
        client.send_message(message)