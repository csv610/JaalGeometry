/****************************************************************************
** Meta object code from reading C++ file 'FontsDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.5.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "FontsDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'FontsDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.5.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JFontsDialog_t {
    QByteArrayData data[7];
    char stringdata0[59];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JFontsDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JFontsDialog_t qt_meta_stringdata_JFontsDialog = {
    {
        QT_MOC_LITERAL(0, 0, 12), // "JFontsDialog"
        QT_MOC_LITERAL(1, 13, 7), // "xRotate"
        QT_MOC_LITERAL(2, 21, 0), // ""
        QT_MOC_LITERAL(3, 22, 7), // "yRotate"
        QT_MOC_LITERAL(4, 30, 7), // "zRotate"
        QT_MOC_LITERAL(5, 38, 11), // "setFontSize"
        QT_MOC_LITERAL(6, 50, 8) // "setColor"

    },
    "JFontsDialog\0xRotate\0\0yRotate\0zRotate\0"
    "setFontSize\0setColor"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JFontsDialog[] = {

// content:
    7,       // revision
    0,       // classname
    0,    0, // classinfo
    5,   14, // methods
    0,    0, // properties
    0,    0, // enums/sets
    0,    0, // constructors
    0,       // flags
    0,       // signalCount

// slots: name, argc, parameters, tag, flags
    1,    0,   39,    2, 0x08 /* Private */,
    3,    0,   40,    2, 0x08 /* Private */,
    4,    0,   41,    2, 0x08 /* Private */,
    5,    0,   42,    2, 0x08 /* Private */,
    6,    0,   43,    2, 0x08 /* Private */,

// slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

    0        // eod
};

void JFontsDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JFontsDialog *_t = static_cast<JFontsDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0:
            _t->xRotate();
            break;
        case 1:
            _t->yRotate();
            break;
        case 2:
            _t->zRotate();
            break;
        case 3:
            _t->setFontSize();
            break;
        case 4:
            _t->setColor();
            break;
        default:
            ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject JFontsDialog::staticMetaObject = {
    {   &QDialog::staticMetaObject, qt_meta_stringdata_JFontsDialog.data,
        qt_meta_data_JFontsDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR
    }
};


const QMetaObject *JFontsDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JFontsDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JFontsDialog.stringdata0))
        return static_cast<void*>(const_cast< JFontsDialog*>(this));
    if (!strcmp(_clname, "Ui::FontsDialog"))
        return static_cast< Ui::FontsDialog*>(const_cast< JFontsDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JFontsDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 5)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 5;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 5)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 5;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
